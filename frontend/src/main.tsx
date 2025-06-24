// import { StrictMode } from 'react';
//
// React StrictMode is a development tool that helps highlight potential problems in an application.
// It activates additional checks and warnings for its descendants. It is ignored in production builds.
//
// Example usage (uncomment for development):
// createRoot(document.getElementById('root')!).render(
//   <StrictMode>
//     <App />
//   </StrictMode>
// )

import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'

createRoot(document.getElementById('root')!).render(
  <App />
)
