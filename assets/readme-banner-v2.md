# README Banner v2 Prompt

Generated asset: `assets/readme-banner-v2.png`

Tool/model: built-in `image_gen` tool

Repo analysis:
- Rust scientific computing crate implementing Minuit-style parameter optimization.
- Primary audience is Rust users doing numerical fitting, scientific computing, and optimization who want a pure-Rust library without C/C++ runtime dependencies.
- Trust surfaces include ROOT Minuit2 parity tests, NIST/real-data examples, CI, cargo-audit, cargo-deny, coverage, and benchmark workflows.
- Visual motifs should come from optimization landscapes, fit curves, covariance ellipses, Hessian/matrix structure, and converging parity traces.

Prompt:

```text
Create a wide 1600x500 GitHub README banner for the open-source Rust scientific computing repository "minuit2-rs".

Style: polished scientific software branding, dark graphite technical background, high contrast, crisp vector-like raster illustration, restrained Rust-orange, cyan, white, and charcoal.

Visual subject: Minuit-style numerical optimization. Show a smooth 3D loss-surface valley, contour rings converging to one minimum, covariance ellipses, fit curves with data points, Hessian/matrix grid motifs, and two aligned orange/cyan optimization traces converging to the same point to suggest ROOT parity.

Text: render exactly "minuit2-rs" as large clean readable text in the left third. No other readable text.

Constraints: no Rust logo, no ROOT or CERN logos, no mascots, no fake badges, no numbers, no slogans, no watermarks.
```

Post-processing:
- Generated source image: `/home/rfrantz/.codex/generated_images/019eb5a2-1499-7ba2-8f25-25684ba2dafd/ig_03393864d525f90a016a2a682babbc8191ae96028aa734f8fe.png`
- Cropped to `1774x554+0+100` and resized to `1600x500` with ImageMagick:
  `magick <source> -crop 1774x554+0+100 +repage -resize 1600x500 assets/readme-banner-v2.png`
