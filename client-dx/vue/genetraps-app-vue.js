const Welcome = {
  template:`
    <v-app>
    <v-btn v-on:click="go()" color="secondary"> Login </v-btn>
    </v-app>`,
    methods: {
            go: function () {
                this.$router.push("/login")
            }
        }
}

const Login = {
  template: `
  <v-app id="login">
    <v-content>
      <v-container fluid fill-height>
        <v-layout column align-center justify-center>
          <v-flex>
            <img src="images/whale.svg" class="float-center" style="maxwidth: 50%; maxheight: 50%;"></img>
          </v-flex>
          <v-flex xs12 sm8 md4>
            <v-card class="elevation-6">
              <v-toolbar dark color="primary">
                <v-toolbar-title>Sign in</v-toolbar-title>
              </v-toolbar>
              <v-card-text>
                <v-form>
                  <v-text-field v-model="login" prepend-icon="person" name="login" label="Login" type="text"></v-text-field>
                  <v-text-field v-model="password" prepend-icon="lock" name="password" label="Password" id="password" type="password"></v-text-field>
                </v-form>
              </v-card-text>
              <v-card-actions>
                <v-spacer></v-spacer>
                <v-btn v-on:click.prevent="getToken" color="primary">Login</v-btn>
              </v-card-actions>
            </v-card>
          </v-flex>
        </v-layout>
      </v-container>
    </v-content>
  </v-app>
  `,
  methods: {
          getToken: function (event) {
              console.log("Getting token")
              console.log(this.login)
              console.log(this.password)

              axios({
                method: 'post', //you can set what request you want to be
                url: 'http://genetraps.intelliseq.pl:8088/oauth/token',
                data: {
                  grant_type: "password",
                  client_id: "web_app",
                  username: "admin",
                  password: "welcome1"
                },
                headers: {
                  'authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
                  'content-type': 'application/x-www-form-urlencoded'
                }
              })
              .then(response => {
                console.log(response)
              })
              .catch(e => {
                console.log(e)
                this.errors.push(e)
              })
          }
      },
      data: function () {
        return {
          login: '',
          password: ''
        }
      },
}

const router = new VueRouter({
  mode: 'history',
  routes: [
    { path: '/', component: Welcome },
    { path: '/login', component: Login }
  ]
})

new Vue({
	router,
  el: '#genetraps-app',
  data: {
    msg: 'Hello World'
  }
})
