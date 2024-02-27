import './App.css'
import './App.scss'
import 'bootstrap/dist/css/bootstrap.min.css'
import Nav from './components/Nav'
import Router from './components/Router'
import Footer from './components/Footer'
import {
  BrowserRouter,
  RouterProvider,
  Routes,
  Route,
  NavLink,
} from "react-router-dom"


function App () {
  return (
    <div className="App">
      <BrowserRouter>
        <Nav />
        <Router />
        <Footer />
      </BrowserRouter>
    </div>
  )
}

export default App
