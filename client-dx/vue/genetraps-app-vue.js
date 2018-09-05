/* event hub */
Vue.prototype.$hub = new Vue();

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
    access_token: ""
  },
  created() {
    console.log("LOG: Vue.App.created()")
    this.$hub.$on('login', (response) => {
      console.log(response)
    });
  },
  mounted() {
    console.log("LOG: Vue.App.mounted()")
    if(this.$cookies.get("refresh_token") == null) {
      
    }
  },
  methods: {

  },
  components: {
    'login-component': Login
  }
})
