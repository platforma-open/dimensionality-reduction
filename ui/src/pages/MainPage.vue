<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAccordionSection, PlAlert, PlBlockPage, PlBtnGhost, PlDropdownRef, PlLogView, PlMaskIcon24, PlNumberField, PlRow, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { ref } from 'vue';

const app = useApp();

// const settingsAreShown = ref(app.model.outputs.UMAPPf === undefined)
const settingsAreShown = ref(true);
const showSettings = () => {
  settingsAreShown.value = true;
};

</script>

<template>
  <PlBlockPage>
    <template #title>Dimensionality Reduction</template>
    <template #append>
      <PlBtnGhost @click.stop="showSettings">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlLogView :log-handle="app.model.outputs.dimReductionLog" label="Log"/>

    <PlSlideModal v-model="settingsAreShown">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
        label="Select dataset"
      />
      <!-- Content hidden until you click ADVANCED SETTINGS -->
      <PlAccordionSection label="ADVANCED SETTINGS">
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
        </PlRow>
        <!-- Add warnings if selected parameters are out of most commonly used bounds -->
        <PlAlert v-if="app.model.args.nPCs < 30" type="warn">
          {{ "Warning: The selected number of PCs is out of the most commonly used range (30-50)" }}
        </PlAlert>
        <PlAlert v-if="app.model.args.nPCs > 50" type="warn">
          {{ "Warning: The selected number of PCs is out of the most commonly used range (30-50)" }}
        </PlAlert>
        <PlAlert v-if="app.model.args.nNeighbors < 10" type="warn">
          {{ "Warning: The selected number of neighbors is out of the most commonly used range (10-30)" }}
        </PlAlert>
        <PlAlert v-if="app.model.args.nNeighbors > 30" type="warn">
          {{ "Warning: The selected number of neighbors is out of the most commonly used range (10-30)" }}
        </PlAlert>
      </PlAccordionSection>
    </PlSlideModal>
  </PlBlockPage>
</template>
