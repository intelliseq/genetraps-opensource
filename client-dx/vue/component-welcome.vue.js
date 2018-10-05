const welcomeComponent = {
  template:`
    <v-app>

    </v-app>`,
    data: function () {
      return {
        waiting_visibility: false
      }
    },
    created: function () {
      logger("DEBUG", "vue.welcome.created")
      //this.$hub.$on('wait', (turnon) => {
        //this.waiting_visibility = turnon
        //console.log("LOG: Vue.Welcome.hub.wait")
      //});
    },
    methods: {
            go: function () {
                this.$router.push("/login")
            }
        }
}
