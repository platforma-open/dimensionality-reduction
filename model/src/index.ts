import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlRef,
} from '@platforma-sdk/model';
import {
  BlockModel,
  isPColumn,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type UiState = {
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
  graphStateUMAPHarmony: GraphMakerState;
  graphStateTSNEHarmony: GraphMakerState;
};

export type BlockArgs = {
  countsRef?: PlRef;
  covariateRefs: PlRef[];
  nPCs: number;
  nNeighbors: number;
  title?: string;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    covariateRefs: [],
    nPCs: 50,
    nNeighbors: 15,
  })

  .withUiState<UiState>({
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
      currentTab: 'settings',
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
      currentTab: null,
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
    graphStateUMAPHarmony: {
      title: 'UMAP (Harmony)',
      template: 'dots',
      currentTab: null,
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
    graphStateTSNEHarmony: {
      title: 'tSNE (Harmony)',
      template: 'dots',
      currentTab: null,
      layersSettings: {
        dots: {
          dotFill: '#99E099',
        },
      },
    },
  })

  .argsValid((ctx) => ctx.args.countsRef !== undefined)

  .output('countsOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/countMatrix' && spec.domain?.['pl7.app/rna-seq/normalized'] === 'false'
    , { includeNativeLabel: false, addLabelAsSuffix: true }),
  )

  .output('metadataOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec) && spec.name === 'pl7.app/metadata'),
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

  .output('UMAPPfPcols', (ctx) => {
    const pCols = ctx.outputs?.resolve('UMAPPf')?.getPColumns();
    if (pCols === undefined)
      return undefined;

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
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

  .output('tSNEPfPcols', (ctx) => {
    const pCols = ctx.outputs?.resolve('tSNEPf')?.getPColumns();
    if (pCols === undefined)
      return undefined;

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('UMAPHarmonyPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('UMAPHarmonyPf')?.getPColumns();
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

  .output('UMAPHarmonyPfPcols', (ctx) => {
    const pCols = ctx.outputs?.resolve('UMAPHarmonyPf')?.getPColumns();
    if (pCols === undefined)
      return undefined;

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('tSNEHarmonyPf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('tSNEHarmonyPf')?.getPColumns();
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

  .output('tSNEHarmonyPfPcols', (ctx) => {
    const pCols = ctx.outputs?.resolve('tSNEHarmonyPf')?.getPColumns();
    if (pCols === undefined)
      return undefined;

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('hasBatchCorrection', (ctx) => {
    return ctx.args.covariateRefs.length > 0;
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Dimensionality Reduction - ${ctx.args.title}`
      : 'Dimensionality Reduction',
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
