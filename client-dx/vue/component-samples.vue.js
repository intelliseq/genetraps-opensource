const samplesComponent = {
  template:`
  <v-layout>
  <v-flex ma-1>
  <v-card>
    <v-list two-line subheader>
      <v-subheader inset>Samples</v-subheader>
            <v-list-tile v-on:click="goToSample($event, sample.id)" v-for="sample in samples" :key="sample.id" avatar>
            <v-list-tile-avatar>
              <v-icon>assignment_ind</v-icon>
            </v-list-tile-avatar>
            <v-list-tile-content>
              <v-list-tile-title>{{ sample.id }}</v-list-tile-title>
              <v-list-tile-sub-title>{{ sample.access }}</v-list-tile-sub-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-btn icon>
                <v-icon color="grey lighten-1">info</v-icon>
              </v-btn>
            </v-list-tile-action>
          </v-list-tile>
        </v-list>
  </v-card>
  </v-flex>
  </v-layout>

    `,
    data: function () {
      return {
        sampletext: "sample-text"
      }
    },
    computed: {
      ...Vuex.mapState('samples', ['samples'])
    },
    methods: {
      createSample: function (event) {
        logger.debug("vue.welcome.createSample")
        store.dispatch('sample/createSample')
      },
      getSamples: function (event) {
        logger.debug("vue.welcome.getSamples")
        store.dispatch('sample/getSamples')
      },
      goToSample: function (event, id){
        router.push({ path: `/sample/${id}`})
      }
    },
    created: function () {
      logger.debug("vue.samples.created")
    }
}
