const store = new Vuex.Store({
    state: {
        waitComponentVisibility: false,
        toolbarVisibility: false,
        buttonColor: 'teal lighten-2',
        buttonColors: {
          'samples':'teal lighten-2',
          'panels':'teal lighten-2',
        },
        waitingText: "waiting",
    },
    mutations: {
      setButtonColors(state, route) {
        if (route) {
          state.buttonColors.samples = "teal lighten-4"
          state.buttonColors.panels = "teal lighten-2"
        }
      },
      setWaitComponentVisibility(state, waitComponentVisibility) {
        state.waitComponentVisibility = waitComponentVisibility
        logger.debug("setWaitComponentVisibility")
      },
      setWaitingText(state, waitingText) {
        state.waitingText = waitingText
      },
      setToolbarVisibility(state, toolbarVisibility) {
        state.toolbarVisibility = toolbarVisibility
        logger.debug("setToolbarVisibility")
      }
    },
    modules: {
    	user: userModule,
      security: securityModule,
      sample: sampleModule,
      samples: samplesModule,
    },
})
