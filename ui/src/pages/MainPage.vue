<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAccordionSection, PlAlert, PlBlockPage, PlDropdownMulti, PlDropdownRef, PlNumberField, PlRow, PlTabs } from '@platforma-sdk/ui-vue';
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

const covariateOptions = computed(() => {
  return app.model.outputs.metadataOptions?.map((v) => ({
    value: v.ref,
    label: v.label,
  })) ?? [];
});

const tabOptions = computed(() => {
  const baseOptions = [
    { label: 'UMAP', value: 'umap' },
    { label: 't-SNE', value: 'tsne' },
  ];

  if (app.model.outputs.hasBatchCorrection) {
    baseOptions.push(
      { label: 'UMAP (Harmony)', value: 'umap-harmony' },
      { label: 't-SNE (Harmony)', value: 'tsne-harmony' },
    );
  }

  return baseOptions;
});

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

/* Function to create default options according to the selected tab */
function createDefaultOptions(
  pcols: PColumnIdAndSpec[] | undefined,
  coord1Name: string,
  coord2Name: string,
): PredefinedGraphOption<'scatterplot-umap'>[] | undefined {
  if (!pcols || pcols.length === 0)
    return undefined;

  const coord1Index = getIndex(coord1Name, pcols);
  const coord2Index = getIndex(coord2Name, pcols);

  if (coord1Index === -1 || coord2Index === -1)
    return undefined;

  return [
    {
      inputName: 'x',
      selectedSource: pcols[coord1Index].spec,
    },
    {
      inputName: 'y',
      selectedSource: pcols[coord2Index].spec,
    },
    /* Grouping/Color is set to sampleId, first axis of coord pcols */
    {
      inputName: 'grouping',
      selectedSource: pcols[0].spec.axesSpec[0],
    },
  ];
}

const defaultOptions = computed((): PredefinedGraphOption<'scatterplot-umap'>[] | undefined => {
  if (data.currentTab === 'umap') {
    return createDefaultOptions(
      app.model.outputs.UMAPPfPcols,
      'pl7.app/rna-seq/umap1',
      'pl7.app/rna-seq/umap2',
    );
  }
  if (data.currentTab === 'tsne') {
    return createDefaultOptions(
      app.model.outputs.tSNEPfPcols,
      'pl7.app/rna-seq/tsne1',
      'pl7.app/rna-seq/tsne2',
    );
  }
  if (data.currentTab === 'umap-harmony') {
    return createDefaultOptions(
      app.model.outputs.UMAPHarmonyPfPcols,
      'pl7.app/rna-seq/umap1',
      'pl7.app/rna-seq/umap2',
    );
  }
  if (data.currentTab === 'tsne-harmony') {
    return createDefaultOptions(
      app.model.outputs.tSNEHarmonyPfPcols,
      'pl7.app/rna-seq/tsne1',
      'pl7.app/rna-seq/tsne2',
    );
  }
  return undefined;
});

/* Modify graph state, pframe and default options based on the selected tab */
const graphState = computed({
  get: () => {
    if (data.currentTab === 'umap') return app.model.ui.graphStateUMAP;
    if (data.currentTab === 'tsne') return app.model.ui.graphStateTSNE;
    if (data.currentTab === 'umap-harmony') return app.model.ui.graphStateUMAPHarmony;
    if (data.currentTab === 'tsne-harmony') return app.model.ui.graphStateTSNEHarmony;
    return app.model.ui.graphStateUMAP;
  },
  set: (value) => {
    if (data.currentTab === 'umap') app.model.ui.graphStateUMAP = value;
    else if (data.currentTab === 'tsne') app.model.ui.graphStateTSNE = value;
    else if (data.currentTab === 'umap-harmony') app.model.ui.graphStateUMAPHarmony = value;
    else if (data.currentTab === 'tsne-harmony') app.model.ui.graphStateTSNEHarmony = value;
  },
});

const pFrame = computed(() => {
  if (data.currentTab === 'umap') return app.model.outputs.UMAPPf;
  if (data.currentTab === 'tsne') return app.model.outputs.tSNEPf;
  if (data.currentTab === 'umap-harmony') return app.model.outputs.UMAPHarmonyPf;
  if (data.currentTab === 'tsne-harmony') return app.model.outputs.tSNEHarmonyPf;
  return app.model.outputs.UMAPPf;
});

/* Use both currentTab and pFrame in :key to force re-render the graph when either args (which changes the pFrame) or the tab changes */

</script>

<template>
  <PlBlockPage>
    <GraphMaker
      :key="`${data.currentTab}-${pFrame}`"
      v-model="graphState"
      chartType="scatterplot-umap"
      :p-frame="pFrame"
      :default-options="defaultOptions"
    >
      <template #titleLineSlot>
        <PlTabs v-model="data.currentTab" :options="tabOptions" :style="{ display: 'flex', justifyContent: 'flex-end' }"/>
      </template>
      <template #settingsSlot>
        <PlDropdownRef
          v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
          :style="{ width: '320px' }"
          label="Select dataset"
          clearable
          required
          @update:model-value="setInput"
        />
        <PlDropdownMulti
          v-model="app.model.args.covariateRefs"
          :options="covariateOptions"
          :style="{ width: '320px' }"
          label="Batch correction covariates (optional)"
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
