const userModule = {
	namespaced: true,
	state: {
		userFirstName: '',
		userLastName: '',
		userName: '',
		userId: '',
		userEmail: ''
  },
	mutations: {
		setUserData (state, data) {
			state.userFirstName = data.FirstName
			state.userLastName = data.LastName
			state.userName = data.Username
			state.userId = data.UserID
			state.userEmail = data.Email
			logger.debug("Set user data for " + state.userEmail)
	  }
	},
  actions: {
    getUser() {
      logger.debug("vue.vuex.user.getUser")
			request({
				waitingText: "Verifying user credentials",
				endpoint: "user",
				callback: function(data) {
					store.commit('user/setUserData', data)
				} // calback:{}
			}) // request()
    }, // getUser()
  } // actions:{}
}
