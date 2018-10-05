const welcomeComponent = {
  template:`
    <v-app>
      <v-btn v-on:click.prevent="createSample">CREATE SAMPLE</v-btn>
    </v-app>`,
    data: function () {
      return {

      }
    },
    methods: {
      createSample: function (event) {
        logger.debug("vue.welcome.createSample")
        store.dispatch('sample/createSample')
      }
    },
    created: function () {
      logger.debug("vue.welcome.created")
    }
}
