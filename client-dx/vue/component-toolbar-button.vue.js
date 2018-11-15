const toolbarButtonComponent = {
  template:`
  
  <v-tooltip right>
  <v-btn depressed dark color="teal lighten-2" slot="activator" block>
    <v-layout row justify-start align-center>
      <v-icon>{{item.icon}}</v-icon>
      <span>&nbsp;&nbsp;{{item.title}}</span>
    </v-layout>
  </v-btn>
  <span>{{item.tooltip}}</span>
  </v-tooltip>

`,
    props: ['item'],
    created: function () {
      logger.debug("vue.toolbarButton.created")
    }
}
