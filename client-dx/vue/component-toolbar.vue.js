const toolbarComponent = {
  template:`

  <v-toolbar absolute color="teal lighten-3" dark scroll-off-screen scroll-target="#scrolling-techniques">
    <v-img src="images/logo-white.svg" height="90%" aspect-ratio="1" contain></v-img>
    <v-toolbar-side-icon></v-toolbar-side-icon>
    <v-toolbar-title>genetraps</v-toolbar-title>
    <v-spacer></v-spacer>
    <span>{{userEmail}}</span>
  </v-toolbar>

`,
    computed: {
      ...Vuex.mapState('user', ['userEmail']),
      ...Vuex.mapState(['toolbarVisibility'])

    },
    //  Vuex.mapState('user', ['userEmail'])
    //}
    created: function () {
      logger.debug("vue.wait.created")
    }
}
