import React from 'react';
import ReactDOM from 'react-dom/client';

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css';

import reportWebVitals from './reportWebVitals';
import { createBrowserRouter, RouterProvider } from 'react-router-dom';


import App from './App';
import Similarity from './Similarity/Similarity';
import Navigation from './Navbar';
import Substructure from './Substructure/Substructure';
import Substances from './Substances/Substances';

const root = ReactDOM.createRoot(document.getElementById('root'));

const router = createBrowserRouter([
  {
    path: "/",
    element: <App />,
  },
  {
    path: "/Similarity",
    element: <Similarity />,
  },
  {
    path: "/Substructure",
    element: <Substructure />,
  },
  {
    path: "/Substances",
    element: <Substances />,
  }
]);



root.render(
  <div>
    <Navigation />
    <RouterProvider router={router} />
  </div>

);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();
