import React from 'react'
import { useRoutes } from 'react-router-dom'
import { Browser } from './Browser'
import Banner from './Banner'
import { Projects } from './Projects'
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
      element: <Projects />,
    },
    {
      path: 'projects',
      element: <Projects />,
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
