import React from 'react'
import { useRoutes } from 'react-router-dom'
import { Tools } from './Tools'
import Banner from './Banner'
import { Projects } from './Projects'
import Browser from './Browser'
import { Contact } from './Contact'
import { Database } from './Database'
import Submit from './Submit'

export default function Router() {
  let routes = useRoutes([
    {
      path: '/',
      element: <Banner />,
    },
    {
      path: '/home',
      element: <Banner />,
    },
    {
      path: 'browse',
      element: <Browser />,
    },
    {
      path: 'tools',
      element: <Tools />,
    },
    {
      path: 'contact',
      element: <Contact />,
    },
    {
      path: 'database',
      element: <Database />,
    },
    {
      path: 'submit',
      element: <Submit />,
      children: [
        {
          path: 'link',
          element: <Submit />,
        },
      ],
    },
  ])
  return routes
}
