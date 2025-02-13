<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { ref } from 'vue';
import { GraphMaker } from '@milaboratories/graph-maker';

const app = useApp();

</script>

<template>
  <PlBlockPage>
    <!-- <template #title>Dimensionality Reduction</template> -->
    <GraphMaker chartType="scatterplot" :p-frame="app.model.outputs.tSNEPf" v-model="app.model.ui.graphState" />
  </PlBlockPage>
</template>