import { model } from '@platforma-open/milaboratories.dimensionality-reduction.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import MainPage from './pages/MainPage.vue';
import UMAP from './pages/UMAP.vue';
import tSNE from './pages/tSNE.vue';

export const sdkPlugin = defineApp(model, (app) => {
  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    showErrorsNotification: true,
    routes: {
      '/': () => MainPage,
      '/umap': () => UMAP,
      '/tsne': () => tSNE,
    },
  };
});

export const useApp = sdkPlugin.useApp;
