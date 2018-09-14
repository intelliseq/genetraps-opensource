const CookiesComponent = {
  template:`
  <v-app>
  <v-layout column align-center justify-center>
    <v-flex text-xs-center transition="scale-transition">
      <img src="images/logo.svg" class="float-center" style="width: 25%;"></img>
    </v-flex>
  </v-layout>
  <v-layout row justify-center>
    <v-dialog v-model="dialog_visibility" persistent max-width="500">
      <v-card>
        <v-card-title class="headline">genomicapt uses cookies</v-card-title>
        <v-card-text>We use cookies for security purposes, to personalize content and to analyse our traffic. We use Google analytics to analyse traffic. You consent to our cookies if you continue to use this website.</v-card-text>
        <v-card-actions>
          <v-layout column>
            <v-spacer></v-spacer>
            <v-layout row>
              <v-checkbox v-model="necessary" disabled color="grey lighten-1" :label="'Necessary'" max-width=100></v-checkbox>
              <v-checkbox v-model="preferences" color="success" :label="'Preferences'" max-width=100></v-checkbox>
              <v-checkbox v-model="statistics" color="success" :label="'Statistics'" max-width=100></v-checkbox>
            </v-layout>
            <v-spacer></v-spacer>
            <v-btn color="green darken-1" flat v-on:click.prevent="setCookies">Agree</v-btn>
          <v-layout>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </v-layout>
  </v-app>`,
  mounted() {
    console.log("LOG: Vue.Cookies.mounted()")
    this.dialog_visibility = true
  },
  data: function () {
    return {
      dialog_visibility: true,
      necessary: true,
      preferences: true,
      statistics: true
    }
  },

    methods: {
            setCookies: function () {
              console.log("LOG: Vue.Cookies.setCookies()")
              console.log("LOG: Vue.Cookies.setCookies() necessary " + this.necessary)
              console.log("LOG: Vue.Cookies.setCookies() prefrences " + this.preferences)
              console.log("LOG: Vue.Cookies.setCookies() statistics " + this.statistics)
              this.$cookies.set("cookies_necessary", true)
              if(this.preferences) {
                console.log("LOG: Vue.Cookies.setCookies() preferences")
                this.$cookies.set("cookies_preferences", true)
              }
              if(this.statistics) {
                console.log("LOG: Vue.Cookies.setCookies() statistics")
                this.$cookies.set("cookies_statistics", true)
                loadGoogleAnalytics()
              }
              this.$router.push("/login")
            }
        }
}
