<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAccordionSection, PlAlert, PlBlockPage, PlDropdownRef, PlNumberField, PlRow, PlTabs } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import type { PColumnIdAndSpec, PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';
import { computed, reactive } from 'vue';

const app = useApp();

const data = reactive({
  currentTab: 'umap',
});

const tabOptions = [
  { label: 'UMAP', value: 'umap' },
  { label: 't-SNE', value: 'tsne' },
];

function setInput(inputRef?: PlRef) {
  app.model.args.countsRef = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.countsOptions?.find((o: { ref: PlRef; label: string }) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

function getIndex(name: string, pcols: PColumnIdAndSpec[]): number {
  return pcols.findIndex((p) => (p.spec.name === name
  ));
}

const defaultOptionsUMAP = computed((): PredefinedGraphOption<'scatterplot-umap'>[] | undefined => {
  const pcols = app.model.outputs.UMAPPfPcols;
  if (!pcols)
    return undefined;

  return [
    {
      inputName: 'x',
      selectedSource: pcols[getIndex('pl7.app/rna-seq/umap1', pcols)].spec,
    },
    {
      inputName: 'y',
      selectedSource: pcols[getIndex('pl7.app/rna-seq/umap2', pcols)].spec,
    },
    /* Grouping/Color is set to sampleId, first axis of coord pcols */
    {
      inputName: 'grouping',
      selectedSource: pcols[0].spec.axesSpec[0],
    },
  ];
});

const defaultOptionsTSNE = computed((): PredefinedGraphOption<'scatterplot-umap'>[] | undefined => {
  const pcols = app.model.outputs.tSNEPfPcols;
  if (!pcols)
    return undefined;

  return [
    {
      inputName: 'x',
      selectedSource: pcols[getIndex('pl7.app/rna-seq/tsne1', pcols)].spec,
    },
    {
      inputName: 'y',
      selectedSource: pcols[getIndex('pl7.app/rna-seq/tsne2', pcols)].spec,
    },
    {
      inputName: 'grouping',
      selectedSource: pcols[0].spec.axesSpec[0],
    },
  ];
});

/* Modify graph state, pframe and default options based on the selected tab */
const graphState = computed({
  get: () => data.currentTab === 'umap' ? app.model.ui.graphStateUMAP : app.model.ui.graphStateTSNE,
  set: (value) => {
    if (data.currentTab === 'umap')
      app.model.ui.graphStateUMAP = value;
    else
      app.model.ui.graphStateTSNE = value;
  },
});

const pFrame = computed(() => data.currentTab === 'umap' ? app.model.outputs.UMAPPf : app.model.outputs.tSNEPf);
const defaultOptions = computed(() => data.currentTab === 'umap' ? defaultOptionsUMAP.value : defaultOptionsTSNE.value);

/* Use both currentTab and pFrame in :key to force re-render the graph when either args (which changes the pFrame) or the tab changes */

</script>

<template>
  <PlBlockPage>
    <PlTabs v-model="data.currentTab" :options="tabOptions" />
    <GraphMaker
      :key="`${data.currentTab}-${pFrame}`"
      v-model="graphState"
      chartType="scatterplot-umap"
      :p-frame="pFrame"
      :default-options="defaultOptions"
    >
      <template #settingsSlot>
        <PlDropdownRef
          v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
          :style="{ width: '320px' }"
          label="Select dataset"
          clearable
          required
          @update:model-value="setInput"
        />
        <!-- Content hidden until you click ADVANCED SETTINGS -->
        <PlAccordionSection :style="{ width: '320px' }" label="ADVANCED SETTINGS">
          <PlRow>
            <PlNumberField
              v-model="app.model.args.nPCs"
              label="N PCs" :min-value="20" :step="1"
            >
              <template #tooltip>
                <div>
                  <strong>Number of Principal Components (PCs)</strong><br/>
                  Controls how many principal components are used for UMAP and t-SNE calculations.<br/><br/>
                  <strong>Recommended ranges:</strong><br/>
                  • 30-50: Optimal for most datasets<br/>
                  • 20-30: For smaller datasets or faster computation<br/>
                  • 50+: For very large datasets with many genes<br/><br/>
                  <strong>Effect:</strong> Higher values capture more variance but increase computation time.
                </div>
              </template>
            </PlNumberField>
            <PlNumberField
              v-model="app.model.args.nNeighbors"
              label="N Neighbors" :min-value="5" :step="1"
            >
              <template #tooltip>
                <div>
                  <strong>Number of Neighbors for UMAP</strong><br/>
                  Controls the balance between local and global structure in UMAP visualization.<br/><br/>
                  <strong>Recommended ranges:</strong><br/>
                  • 10-30: Optimal for most datasets<br/>
                  • 5-10: Emphasizes local structure (more clusters)<br/>
                  • 30+: Emphasizes global structure (fewer clusters)<br/><br/>
                  <strong>Effect:</strong> Lower values preserve local structure, higher values preserve global relationships.
                </div>
              </template>
            </PlNumberField>
          </PlRow>
          <!-- Add warnings if selected parameters are out of most commonly used bounds -->
          <PlAlert v-if="app.model.args.nPCs > 20 && app.model.args.nPCs < 30" type="warn">
            <template #title>Suboptimal PC Count</template>
            The selected number of PCs ({{ app.model.args.nPCs }}) is below the recommended range (30-50).
            This may result in loss of important biological variation. Consider increasing to 30-50 for better results.
          </PlAlert>
          <PlAlert v-if="app.model.args.nPCs > 50" type="warn">
            <template #title>High PC Count</template>
            The selected number of PCs ({{ app.model.args.nPCs }}) is above the recommended range (30-50).
            This will increase computation time without significant improvement in results for most datasets.
          </PlAlert>
          <PlAlert v-if="app.model.args.nNeighbors > 5 && app.model.args.nNeighbors < 10" type="warn">
            <template #title>Low Neighbor Count</template>
            The selected number of neighbors ({{ app.model.args.nNeighbors }}) is below the recommended range (10-30).
            This may result in fragmented clusters and loss of global structure.
          </PlAlert>
          <PlAlert v-if="app.model.args.nNeighbors > 30" type="warn">
            <template #title>High Neighbor Count</template>
            The selected number of neighbors ({{ app.model.args.nNeighbors }}) is above the recommended range (10-30).
            This may over-smooth the visualization and obscure fine-grained cell type distinctions.
          </PlAlert>
        </PlAccordionSection>
      </template>
    </GraphMaker>
  </PlBlockPage>
</template>
