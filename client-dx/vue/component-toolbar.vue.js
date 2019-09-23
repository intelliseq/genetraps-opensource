const toolbarComponent = {
  template:`
  <v-toolbar color="teal lighten-3" dark height="50px" max-height="50px" style="padding: 0!important">

    <v-flex xs12 sm2 md2 hidden-xs-only>
      <v-img v-on:click="setToolbarVisibilityToFalse()" src="images/logo-white.svg" max-height="40px" aspect-ratio="1" contain></v-img>
    </v-flex>
    <!--store.commit("setToolbarVisibility", true)-->
    <!--window.location.href = 'http://localhost:8080'-->

    <v-btn v-on:click.prevent="showSamples" :color="buttonColors.samples" depressed>
      <v-flex hidden-sm-and-down><span>SAMPLES</span></v-flex>
      <v-flex hidden-sm-and-down><v-icon right dark>folder_shared</v-icon></v-flex>
      <v-flex hidden-md-and-up><v-icon dark>folder_shared</v-icon></v-flex>
    </v-btn>

    <v-btn :color="buttonColors.panels" depressed>
      <v-flex hidden-sm-and-down><span>PANELS</span></v-flex>
      <v-flex hidden-sm-and-down><v-icon dark right>assignment</v-icon></v-flex>
      <v-flex hidden-md-and-up><v-icon dark>assignment</v-icon></v-flex>
    </v-btn>

    <v-spacer></v-spacer>

    <v-btn color="teal lighten-2" depressed>
      <v-flex hidden-sm-and-down><span>{{userEmail}}</span></v-flex>
      <v-flex hidden-sm-and-down><v-icon right dark>arrow_drop_down_circle</v-icon></v-flex>
      <v-flex hidden-md-and-up><v-icon dark>arrow_drop_down_circle</v-icon></v-flex>
    </v-btn>

  </v-toolbar>
  `,
    computed: {
      ...Vuex.mapState('user', ['userEmail']),
      ...Vuex.mapState(['buttonColors'])

    },
    methods: {
      showSamples: function (event) {
        logger.debug("vue.toolbar.showSamples")
        store.dispatch('samples/showSamples')
      },
      setToolbarVisibilityToFalse(event) {
        store.commit('setToolbarVisibility', false)
        router.push('/')
      }
    },
    created: function () {
      logger.debug("vue.wait.created")
    }
}
