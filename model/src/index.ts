import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PFrameHandle,
  PlRef } from '@platforma-sdk/model';
import {
  BlockModel,
  createPlDataTable,
  isPColumn,
  isPColumnSpec,
  PlDataTableState,
} from '@platforma-sdk/model';

export type UiState = {
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
};

// export type Formula = {
//   // we put formula label in the arg as it will be used
//   // in the annotations to re-use in the downstream blocks
//   label: string;
//   covariateRefs: PlRef[];
//   contrastFactor?: PlRef;
//   denominator?: String;
//   numerator?: String;
// };

export type BlockArgs = {
  countsRef?: PlRef;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
  })

  .withUiState<UiState>({
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
    },
  })

  .output('countsOptions', (ctx) =>
    // I've added these "||" for backward compatibility (As I see, the shape of PColum was changed)
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/countMatrix' && spec.domain?.['pl7.app/rna-seq/normalized'] === 'false'
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('UMAPPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('UMAPPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    // enriching with upstream data
    const upstream = ctx.resultPool
      .getData()
      .entries.map((v) => v.obj)
      .filter(isPColumn)
      .filter((column) => column.id.includes('metadata'));

    return ctx.createPFrame([...pCols, ...upstream]);
  })

  .output('tSNEPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('tSNEPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    // enriching with upstream data
    const upstream = ctx.resultPool
      .getData()
      .entries.map((v) => v.obj)
      .filter(isPColumn)
      .filter((column) => column.id.includes('metadata'));

    return ctx.createPFrame([...pCols, ...upstream]);
  })

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/tsne', label: 'tSNE' },
  ]))

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
