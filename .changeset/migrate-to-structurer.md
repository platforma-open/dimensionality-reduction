---
'@platforma-open/milaboratories.dimensionality-reduction.model': patch
'@platforma-open/milaboratories.dimensionality-reduction.ui': patch
'@platforma-open/milaboratories.dimensionality-reduction.workflow': patch
'@platforma-open/milaboratories.dimensionality-reduction.software': patch
'@platforma-open/milaboratories.dimensionality-reduction': patch
---

Migrate block onto the structurer (block-tools 2.12.13) — full SDK upgrade: model/ui-vue 1.81.1, workflow-tengo 6.8.2, tengo-builder 4.0.22, package-builder 3.14.2, test 1.81.2. Adopts the canonical tool-managed layout (oxlint/oxfmt, tsconfig, turbo, CI workflows, managed package.json + catalog) and the slim facade for the root block package. Author-code fixes for the SDK majors: explicit type argument on the `isPColumn` filter feeding `createPFrame`, removal of the retired `@platforma-sdk/ui-vue/styles` import, and the model export renamed `model` -> `platforma` for the facade.
