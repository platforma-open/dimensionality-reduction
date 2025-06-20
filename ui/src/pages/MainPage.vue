<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAccordionSection, PlAlert, PlBlockPage, PlDropdownRef, PlNumberField, PlRow } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import type { GraphMakerProps } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import { plRefsEqual, type PlRef } from '@platforma-sdk/model';

const app = useApp();

function setInput(inputRef?: PlRef) {
  app.model.args.countsRef = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.countsOptions?.find((o: { ref: PlRef; label: string }) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

const defaultOptions: GraphMakerProps['defaultOptions'] = [
  {
    inputName: 'x',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/umap1',
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
      name: 'pl7.app/rna-seq/umap2',
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
    <GraphMaker
      v-model="app.model.ui.graphStateUMAP"
      chartType="scatterplot-umap"
      :p-frame="app.model.outputs.UMAPPf"
      :default-options="defaultOptions"
    >
      <template #settingsSlot>
        <PlDropdownRef
          v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
          :style="{ width: '320px' }"
          label="Select dataset"
          clearable @update:model-value="setInput"
        />
        <!-- Content hidden until you click ADVANCED SETTINGS -->
        <PlAccordionSection :style="{ width: '320px' }" label="ADVANCED SETTINGS">
          <PlRow>
            <PlNumberField
              v-model="app.model.args.nPCs"
              label="N PCs" :minValue="20" :step="1"
            >
              <template #tooltip>
                Select number of Principal Components (PCs) to use.
              </template>
            </PlNumberField>
            <PlNumberField
              v-model="app.model.args.nNeighbors"
              label="N Neighbors" :minValue="5" :step="1"
            >
              <template #tooltip>
                Select number of neighbors for the UMAP algorithm.
              </template>
            </PlNumberField>
          </PlRow >
          <!-- Add warnings if selected parameters are out of most commonly used bounds -->
          <PlAlert v-if="app.model.args.nPCs > 20 && app.model.args.nPCs < 30" type="warn">
            {{ "Warning: The selected number of PCs is out of the most commonly used range (30-50)" }}
          </PlAlert>
          <PlAlert v-if="app.model.args.nPCs > 50" type="warn">
            {{ "Warning: The selected number of PCs is out of the most commonly used range (30-50)" }}
          </PlAlert>
          <PlAlert v-if="app.model.args.nNeighbors > 5 && app.model.args.nNeighbors < 10" type="warn">
            {{ "Warning: The selected number of neighbors is out of the most commonly used range (10-30)" }}
          </PlAlert>
          <PlAlert v-if="app.model.args.nNeighbors > 30" type="warn">
            {{ "Warning: The selected number of neighbors is out of the most commonly used range (10-30)" }}
          </PlAlert>
        </PlAccordionSection>
      </template>
    </GraphMaker>
  </PlBlockPage>
</template>
