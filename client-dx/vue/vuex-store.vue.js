const store = new Vuex.Store({
    state: {
        waitingVisibility: false,
        toolbarVisibility: true,
        userEmail: "piechota@gmail.com",
        waitingText: "waiting"
    },
    mutations: {
      setWaitingVisibility(state, waitingVisibility) {
        state.waitingVisibility = waitingVisibility
        logger("DEBUG", "setWaitingVisibility")
      },
      setWaitingText(state, waitingText) {
        state.waitingText = waitingText
      }
    },
    modules: {
    	user: userModule,
      security: securityModule
    },
})
