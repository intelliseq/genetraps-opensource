const toolbarComponent = {
  template:`
  <v-toolbar color="teal lighten-3" dark height="50px" max-height="50px" style="padding: 0!important">

    <v-flex xs12 sm2 md2 hidden-xs-only>
      <v-img v-on:click="window.location.href = 'http://localhost:8080'" src="images/logo-white.svg" max-height="40px" aspect-ratio="1" contain></v-img>
    </v-flex>

    <v-btn :color="buttonColors.samples" depressed>
      <v-flex hidden-sm-and-down><span>SAMPLES</span></v-flex>
      <v-icon right dark>folder_shared</v-icon>
    </v-btn>

    <v-btn :color="buttonColors.panels" depressed>
      <v-flex hidden-sm-and-down><span>PANELS</span></v-flex>
      <v-icon dark right>assignment</v-icon>
    </v-btn>

    <v-spacer></v-spacer>

    <v-btn color="teal lighten-2" depressed>
      <v-flex hidden-sm-and-down><span>{{userEmail}}</span></v-flex>
      <v-icon right dark>arrow_drop_down_circle</v-icon>
    </v-btn>

  </v-toolbar>
  `,
    computed: {
      ...Vuex.mapState('user', ['userEmail']),
      ...Vuex.mapState(['buttonColors'])

    },
    //  Vuex.mapState('user', ['userEmail'])
    //}
    created: function () {
      logger.debug("vue.wait.created")
    }
}
