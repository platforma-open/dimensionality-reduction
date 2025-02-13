<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { ref } from 'vue';
import { GraphMaker } from '@milaboratories/graph-maker';

const app = useApp();

const settingsAreShown = ref(app.model.outputs.UMAPPf === undefined)
const showSettings = () => { settingsAreShown.value = true }

</script>

<template>
  <PlBlockPage>

    <!-- <template #title>Dimensionality Reduction</template> -->
    <template #append>
      <PlBtnGhost @click.stop="showSettings">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlSlideModal v-model="settingsAreShown">
      <template #title>Settings</template>
      <PlDropdownRef v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
        label="Select dataset" />
    </PlSlideModal>

  </PlBlockPage>
</template>
