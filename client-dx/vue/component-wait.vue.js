const waitComponent = {
  template:`
    <v-layout row justify-center>
      <v-dialog v-model="waitComponentVisibility" hide-overlay persistent width="300">
        <v-card color="primary" dark>
          <v-card-text>
            {{waitingText}}
            <v-progress-linear indeterminate color="white" class="mb-0"></v-progress-linear>
          </v-card-text>
        </v-card>
      </v-dialog>
    </v-layout>`,
    computed: Vuex.mapState(['waitComponentVisibility','waitingText']),
    created: function () {
      logger.debug("vue.wait.created")
    }
}
