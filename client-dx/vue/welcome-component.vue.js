const Wait = {
  template:`
    <v-layout row justify-center>
      <v-dialog v-model="waiting_visibility" hide-overlay persistent width="300">
        <v-card color="primary" dark>
          <v-card-text>
            {{waiting_text}}
            <v-progress-linear indeterminate color="white" class="mb-0"></v-progress-linear>
          </v-card-text>
        </v-card>
      </v-dialog>
    </v-layout>`,
    data: function () {
      return {
        waiting_text: "Signing in",
        waiting_visibility: false
      }
    },
    created: function () {
      console.log("LOG: Vue.Welcome.created()")
      this.$hub.$on('start_waiting', (text) => {
        this.waiting_text = text
        this.waiting_visibility = true
        console.log("LOG: Vue.Welcome.hub.wait")
      });
      this.$hub.$off('stop_waiting', () => {
        this.waiting_visibility = false
        console.log("LOG: Vue.Welcome.hub.wait")
      });
    }
}

const Welcome = {
  template:`
    <v-app>
    <wait></wait>
    </v-app>`,
    data: function () {
      return {
        waiting_visibility: false
      }
    },
    created: function () {
      console.log("LOG: Vue.Welcome.created()")
      this.$hub.$on('wait', (turnon) => {
        this.waiting_visibility = turnon
        console.log("LOG: Vue.Welcome.hub.wait")
      });
    },
    components: {
      "wait": Wait
    },
    methods: {
            go: function () {
                this.$router.push("/login")
            }
        }
}
