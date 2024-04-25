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
import { ConfigProvider, theme } from 'antd'


function App() {
  return (
    <ConfigProvider theme={{
      cssVar: true,
      hashed: false,
      algorithm: theme.darkAlgorithm,
      token: {
        // Seed Token，影响范围大
        colorPrimary: '#00b96b',
        borderRadius: 2,

        // 派生变量，影响范围小
        colorBgContainer: '#f6ffed',
      },
    }}>
      <div className="App">
        <BrowserRouter>
          <Nav />
          <Router />
          <Footer />
        </BrowserRouter>
      </div>
    </ConfigProvider>
  )
}

export default App
