const samplesComponent = {
  template:`

  <div class="text-xs-center">
    <v-btn v-on:click.prevent="createSample">CREATE SAMPLE</v-btn>
  </div>

    `,
    data: function () {
      return {

      }
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
    },
    created: function () {
      logger.debug("vue.welcome.created")
    }
}
