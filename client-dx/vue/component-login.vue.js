const loginComponent = {
  template: `
  <v-app>
  <v-alert :value="error" type="error" dismissible transition="scale-transition" v-model="alert_visibility">
    {{error_message}}
    <br/>
    {{error_description}}
  </v-alert>

        <v-layout column align-center justify-center>
          <v-flex text-xs-center transition="scale-transition">
            <img src="images/logo.svg" class="float-center" style="width: 25%;"></img>
          </v-flex>
          <v-flex>
            <v-dialog v-model="dialog_visibility" persistent max-width="500px">
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
            </v-dialog>
          </v-flex>
        </v-layout>
        </v-app>
  `,
  methods: {
          getToken: function (event) {
              logger("DEBUG", "vue.login.getToken")
              var credentials = {login: this.login, password: this.password}
              store.dispatch('security/loginWithCredentials', credentials)


          }
      },
      data: function () {
        return {
          login: '',
          password: '',
          dialog_visibility: false,
          error: false,
          error_message: '',
          error_description: '',
          alert_visibility: false
        }
      },
      watch: {
        alert_visibility (val) {
          console.log("LOG: Vue.Login.watch().alert_visibility("+val+")")
          if (!val) {
            this.dialog_visibility = true
          }
        }
      },
      mounted() {
        logger("DEBUG", "vue.login.mounted")
        this.dialog_visibility = true
      }
}