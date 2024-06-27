import './App.css'
import './App.scss'
import 'bootstrap/dist/css/bootstrap.min.css'
import Nav from './components/Nav'
import Router from './components/Router'
import {
  BrowserRouter,
  RouterProvider,
  Routes,
  Route,
  NavLink,
} from "react-router-dom"
import { ConfigProvider, theme } from 'antd'


function App () {
  return (
    <ConfigProvider theme={{
      cssVar: true,
      hashed: false,
      token: {
        // Seed Token，影响范围大
        colorPrimary: '#4ab311',
        borderRadius: 2,
        colorTextBase: '#121212',
        colorTextLabel: '#121212',
        colorBgBase: '#fff',
        colorBgElevated: '#f8fff3',
        fontSize: 16,
        // 派生变量，影响范围小
      },
      components: {
        Menu: {
          colorBgContainer: '#f8fff3',
        },
        Layout: {
          siderBg: '#f8fff3',
          triggerBg: '#4ab322',
          headerBg: '#006600',
        }
      }
    }}>
      <div className="App">
        <BrowserRouter>
          <Nav />
          <Router />
        </BrowserRouter>
      </div>
    </ConfigProvider>
  )
}

export default App
