<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { GraphMaker, GraphMakerProps } from '@milaboratories/graph-maker';

const app = useApp();

const defaultOptions: GraphMakerProps['defaultOptions'] = [
  {
    inputName: 'x',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/tsne1',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/sampleId',
          type: 'String',
        },
        {
          name: 'pl7.app/cellId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'y',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/tsne2',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/sampleId',
          type: 'String',
        },
        {
          name: 'pl7.app/cellId',
          type: 'String',
        },
      ],
    },
  },
];

</script>

<template>
  <PlBlockPage>
    <template #title>Dimensionality Reduction</template>
    <GraphMaker chartType="scatterplot-umap" :p-frame="app.model.outputs.tSNEPf" v-model="app.model.ui.graphStateTSNE" :default-options="defaultOptions" />
  </PlBlockPage>
</template>