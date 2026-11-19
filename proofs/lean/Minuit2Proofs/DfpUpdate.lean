module

public import Mathlib.LinearAlgebra.Matrix.PosDef
public import Mathlib.Tactic

@[expose] public section

open Matrix

namespace Minuit2Proofs

noncomputable section

variable {n : ℕ}

/-- A scaled symmetric outer product. -/
def rankOne (c : ℝ) (u : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  c • vecMulVec u u

lemma rankOne_mulVec (c : ℝ) (u x : Fin n → ℝ) :
    rankOne c u *ᵥ x = (c * (u ⬝ᵥ x)) • u := by
  ext i
  simp only [rankOne, mulVec, dotProduct, Matrix.smul_apply, vecMulVec_apply,
    smul_eq_mul, Pi.smul_apply, mul_assoc, ← Finset.mul_sum]
  ring

lemma rankOne_quad (c : ℝ) (u x : Fin n → ℝ) :
    x ⬝ᵥ (rankOne c u *ᵥ x) = c * (x ⬝ᵥ u) ^ 2 := by
  rw [rankOne_mulVec, dotProduct_smul, dotProduct_comm u x]
  simp only [smul_eq_mul]
  ring

lemma rankOne_hermitian (c : ℝ) (u : Fin n → ℝ) :
    (rankOne c u).IsHermitian := by
  apply Matrix.IsHermitian.ext
  intro i j
  simp [rankOne, vecMulVec, mul_comm]

/-- Copies src/migrad/builder.rs:431-475; see Nocedal and Wright, chapter 6. -/
def dfpUpdate (V : Matrix (Fin n) (Fin n) ℝ) (dx dg : Fin n → ℝ) :=
  let delgam := dx ⬝ᵥ dg
  let vg := V *ᵥ dg
  let gvg := dg ⬝ᵥ vg
  let flnu := delgam⁻¹ • dx - gvg⁻¹ • vg
  V + rankOne delgam⁻¹ dx - rankOne gvg⁻¹ vg +
    if gvg < delgam then rankOne gvg flnu else 0

lemma curvature_pos {V : Matrix (Fin n) (Fin n) ℝ} (hV : V.PosDef)
    {dx dg : Fin n → ℝ} (hc : dx ⬝ᵥ dg ≠ 0) : 0 < dg ⬝ᵥ (V *ᵥ dg) := by
  have hg : dg ≠ 0 := by
    intro hg
    simp [hg] at hc
  simpa using hV.dotProduct_mulVec_pos hg

lemma dfpUpdate_secant_of_ne_zero {V : Matrix (Fin n) (Fin n) ℝ}
    (hV : V.PosDef) {dx dg : Fin n → ℝ} (hc : dx ⬝ᵥ dg ≠ 0) :
    dfpUpdate V dx dg *ᵥ dg = dx := by
  have hb := (curvature_pos hV hc).ne'
  have hf : ((dx ⬝ᵥ dg)⁻¹ • dx - (dg ⬝ᵥ (V *ᵥ dg))⁻¹ • (V *ᵥ dg)) ⬝ᵥ dg = 0 := by
    simp [sub_dotProduct, smul_dotProduct, dotProduct_comm (V *ᵥ dg) dg,
      inv_mul_cancel₀ hc, inv_mul_cancel₀ hb]
  dsimp only [dfpUpdate]
  split_ifs <;>
    simp [add_mulVec, sub_mulVec, rankOne_mulVec, hf,
      dotProduct_comm (V *ᵥ dg) dg, inv_mul_cancel₀ hc, inv_mul_cancel₀ hb]

/-- The updated matrix satisfies the secant equation. -/
theorem dfpUpdate_secant {V : Matrix (Fin n) (Fin n) ℝ} (hV : V.PosDef)
    {dx dg : Fin n → ℝ} (hc : 0 < dx ⬝ᵥ dg) :
    dfpUpdate V dx dg *ᵥ dg = dx :=
  dfpUpdate_secant_of_ne_zero hV hc.ne'

lemma dfpUpdate_quad (V : Matrix (Fin n) (Fin n) ℝ) (dx dg x : Fin n → ℝ) :
    x ⬝ᵥ (dfpUpdate V dx dg *ᵥ x) =
      x ⬝ᵥ (V *ᵥ x) + (dx ⬝ᵥ dg)⁻¹ * (x ⬝ᵥ dx) ^ 2 -
        (dg ⬝ᵥ (V *ᵥ dg))⁻¹ * (x ⬝ᵥ (V *ᵥ dg)) ^ 2 +
        if dg ⬝ᵥ (V *ᵥ dg) < dx ⬝ᵥ dg then
          (dg ⬝ᵥ (V *ᵥ dg)) *
            (x ⬝ᵥ ((dx ⬝ᵥ dg)⁻¹ • dx -
              (dg ⬝ᵥ (V *ᵥ dg))⁻¹ • (V *ᵥ dg))) ^ 2
        else 0 := by
  dsimp only [dfpUpdate]
  split_ifs <;> simp [add_mulVec, sub_mulVec, dotProduct_add, dotProduct_sub, rankOne_quad]

/-- Positive curvature preserves positive definiteness in both update branches. -/
theorem dfpUpdate_posDef {V : Matrix (Fin n) (Fin n) ℝ} (hV : V.PosDef)
    {dx dg : Fin n → ℝ} (hc : 0 < dx ⬝ᵥ dg) : (dfpUpdate V dx dg).PosDef := by
  have hb := curvature_pos hV hc.ne'
  have hsym (x y : Fin n → ℝ) : x ⬝ᵥ (V *ᵥ y) = y ⬝ᵥ (V *ᵥ x) := by
    simpa using hV.isHermitian.star_dotProduct_mulVec_comm x y
  apply Matrix.PosDef.of_dotProduct_mulVec_pos
  · dsimp only [dfpUpdate]
    split_ifs
    · exact ((hV.isHermitian.add (rankOne_hermitian _ _)).sub
        (rankOne_hermitian _ _)).add (rankOne_hermitian _ _)
    · simpa using (hV.isHermitian.add (rankOne_hermitian _ _)).sub
        (rankOne_hermitian _ _)
  · intro x hx
    simp only [star_trivial]
    let b := dg ⬝ᵥ (V *ᵥ dg)
    let a := x ⬝ᵥ (V *ᵥ dg)
    let y := x - (a / b) • dg
    have hb' : b ≠ 0 := hb.ne'
    have hy : y ⬝ᵥ (V *ᵥ y) = x ⬝ᵥ (V *ᵥ x) - b⁻¹ * a ^ 2 := by
      dsimp [y]
      simp only [mulVec_sub, mulVec_smul, sub_dotProduct, dotProduct_sub,
        smul_dotProduct, dotProduct_smul, smul_eq_mul, ← hsym x dg]
      dsimp [a, b] at *
      field_simp
      ring
    have hcorr : 0 ≤ if b < dx ⬝ᵥ dg then
        b * (x ⬝ᵥ ((dx ⬝ᵥ dg)⁻¹ • dx - b⁻¹ • (V *ᵥ dg))) ^ 2 else 0 := by
      split_ifs
      · exact mul_nonneg hb.le (sq_nonneg _)
      · exact le_rfl
    rw [dfpUpdate_quad]
    change 0 < x ⬝ᵥ (V *ᵥ x) + (dx ⬝ᵥ dg)⁻¹ * (x ⬝ᵥ dx) ^ 2 - b⁻¹ * a ^ 2 + _
    have hstep : 0 ≤ (dx ⬝ᵥ dg)⁻¹ * (x ⬝ᵥ dx) ^ 2 :=
      mul_nonneg (inv_nonneg.mpr hc.le) (sq_nonneg _)
    by_cases hy0 : y = 0
    · have hxeq : x = (a / b) • dg := sub_eq_zero.mp hy0
      have ha : a / b ≠ 0 := by
        intro ha
        apply hx
        rw [hxeq, ha, zero_smul]
      have hxdx : x ⬝ᵥ dx ≠ 0 := by
        rw [hxeq, smul_dotProduct, dotProduct_comm dg dx]
        exact mul_ne_zero ha hc.ne'
      have hpos := mul_pos (inv_pos.mpr hc) (sq_pos_of_ne_zero hxdx)
      have hyzero : y ⬝ᵥ (V *ᵥ y) = 0 := by simp [hy0]
      linarith
    · have hpos : 0 < y ⬝ᵥ (V *ᵥ y) := by
        simpa using hV.dotProduct_mulVec_pos hy0
      linarith

/-- Negative curvature gives a negative quadratic value and excludes positive semidefiniteness. -/
theorem dfpUpdate_neg_delgam_not_posDef {V : Matrix (Fin n) (Fin n) ℝ}
    (hV : V.PosDef) {dx dg : Fin n → ℝ} (hc : dx ⬝ᵥ dg < 0) :
    dg ⬝ᵥ (dfpUpdate V dx dg *ᵥ dg) = dx ⬝ᵥ dg ∧
      ¬ (dfpUpdate V dx dg).PosSemidef := by
  have hs := dfpUpdate_secant_of_ne_zero hV hc.ne
  have hq : dg ⬝ᵥ (dfpUpdate V dx dg *ᵥ dg) = dx ⬝ᵥ dg := by
    rw [hs, dotProduct_comm]
  refine ⟨hq, ?_⟩
  intro hp
  have hn := hp.dotProduct_mulVec_nonneg dg
  simp only [star_trivial] at hn
  rw [hq] at hn
  exact (not_le_of_gt hc) hn

end

end Minuit2Proofs
