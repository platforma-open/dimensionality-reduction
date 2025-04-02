import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PFrameHandle,
  PlRef } from '@platforma-sdk/model';
import {
  BlockModel,
  isPColumn,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type UiState = {
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
};

export type BlockArgs = {
  countsRef?: PlRef;
  nPCs: number;
  nNeighbors: number;
  title?: string;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    nPCs: 50,
    nNeighbors: 15,
  })

  .withUiState<UiState>({
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
  })

  .output('countsOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/countMatrix' && spec.domain?.['pl7.app/rna-seq/normalized'] === 'false'
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('dimReductionLog', (wf) => {
    const logHandle = wf.outputs?.resolve('dimReductionLog');
    if (logHandle === undefined) {
      return undefined;
    }
    return logHandle.getLogHandle();
  })

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

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/tsne', label: 'tSNE' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Dimensionality Reduction - ${ctx.args.title}`
      : 'Dimensionality Reduction',
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
