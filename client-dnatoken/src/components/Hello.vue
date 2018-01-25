<template>
  <div id="hello">
    <h1>{{ msg }}</h1>
    <button v-on:click="getData()">Add 1</button>
    <h2>{{ data }}</h2>
  </div>
</template>

<script>
import axios from 'axios';

export default {
  name: 'hello',
  data: function() {
    return {
      msg: 'Welcome to Explorare Client',
      data: 'Empty'
    }
  },
  methods: {
    getData() {
      axios.get(`https://swapi.co/api/people/1/?format=json`)
      .then(response => {
        // JSON responses are automatically parsed.
        this.data = response.data
      })
      .catch(e => {
        this.errors.push(e)
      })
    }
  }

}
</script>

<style scoped>
#hello {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  margin-top: 60px;
}

h1, h2 {
  font-weight: normal;
}

ul {
  list-style-type: none;
  padding: 0;
}

li {
  display: inline-block;
  margin: 0 10px;
}

a {
  color: #42b983;
}
</style>
