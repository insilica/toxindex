import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

// https://vite.dev/config/
export default defineConfig({
  base: './', // <-- Ensures relative asset paths for GCS
  plugins: [react(), tailwindcss()],
  server: {
    proxy: {
      '/api': 'http://localhost:6513',
      '/socket.io': {
        target: 'http://localhost:6513',
        ws: true,
      },
    },
  }
})  
