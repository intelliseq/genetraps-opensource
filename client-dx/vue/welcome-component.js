const Welcome = {
  template:`
    <v-app>
    <!--<v-btn v-on:click="go()" color="secondary"> Login </v-btn>-->
    </v-app>`,
    methods: {
            go: function () {
                this.$router.push("/login")
            }
        }
}
