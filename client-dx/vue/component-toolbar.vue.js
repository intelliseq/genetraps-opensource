const toolbarComponent = {
  template:`
  <v-card
      color="grey lighten-4"
      flat
      height="200px"
      tile
    >
  <v-toolbar absolute color="teal lighten-3" dark scroll-off-screen scroll-target="#scrolling-techniques">
    <v-toolbar-side-icon></v-toolbar-side-icon>
    <v-toolbar-title>genetraps</v-toolbar-title>
    <v-spacer></v-spacer>
    <span>{{toolbarVisibility}}</span><span>=====</span><span>{{userEmail}}</span>
  </v-toolbar>
</v-card>
`,
    computed: {
      ...Vuex.mapState('user', ['userEmail']),
      ...Vuex.mapState(['toolbarVisibility'])

    },
    //  Vuex.mapState('user', ['userEmail'])
    //}
    created: function () {
      logger("DEBUG", "vue.wait.created")
    }
}
