const sampleComponent = {
  template:`
  <v-layout>
  <v-flex ma-1>
  <v-card>
    <v-list two-line subheader>
      <v-subheader inset>Samples</v-subheader>
        <v-list-tile v-for="(value, key) in sample.properties" :key="key" avatar @click="">
            <v-list-tile-avatar>
              <v-icon>label</v-icon>
            </v-list-tile-avatar>
            <v-list-tile-content>
              <v-list-tile-title>{{ key }}</v-list-tile-title>
              <v-list-tile-sub-title>{{ value }}</v-list-tile-sub-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-btn icon>
                <v-icon color="grey lighten-1">info</v-icon>
              </v-btn>
            </v-list-tile-action>
          </v-list-tile>
          <v-list-tile>
          <v-list-tile-avatar>
            <v-icon color="green darken-2">add_circle</v-icon>
          </v-list-tile-avatar>
          <v-text-field @keyup.enter="add" v-model="newKey" label="Property Name"></v-text-field></td>
          <v-text-field @keyup.enter="add" v-model="newValue" label="Property Value"></v-text-field>
          </v-list-tile>
        </v-list>
  </v-card>
  </v-flex>
  </v-layout>

    `,
    data: function () {
      return {
        newKey: undefined,
        newValue: undefined,
      }
    },
    computed: {
      ...Vuex.mapState('sample', ['sample'])
    },
    methods: {
      add(){
    	   Vue.set(this.sample.properties, this.newKey, this.newValue)
         this.newName = undefined;
         this.newNumber = undefined;
         console.log(this.sample)
      },
      putProperties: function (event) {
        //logger.debug("vue.welcome.createSample")
        //store.dispatch('sample/createSample')
      },

    },
    created: function () {
      logger.debug("vue.sample.created")
    }
}
