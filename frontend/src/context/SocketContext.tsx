import React, { createContext, useContext, useEffect, useRef, useState } from 'react';
import { io, Socket } from 'socket.io-client';

interface SocketContextType {
  socket: Socket | null;
  isConnected: boolean;
}

const SocketContext = createContext<SocketContextType>({
  socket: null,
  isConnected: false,
});

export const useSocket = () => {
  const context = useContext(SocketContext);
  if (!context) {
    throw new Error('useSocket must be used within a SocketProvider');
  }
  return context;
};

interface SocketProviderProps {
  children: React.ReactNode;
}

export const SocketProvider: React.FC<SocketProviderProps> = ({ children }) => {
  const socketRef = useRef<Socket | null>(null);
  const [isConnected, setIsConnected] = useState(false);

  useEffect(() => {
    // Create socket only once
    if (!socketRef.current) {
      socketRef.current = io('/', {
        withCredentials: true,
        transports: ['polling'],
        transportOptions: {
          polling: {
            timeout: 30000 // 30 seconds
          }
        },
        // Add connection options to reduce polling frequency
        reconnectionAttempts: 5,
        reconnectionDelay: 1000,
        reconnectionDelayMax: 5000
      });

      // Add connection event handlers
      socketRef.current.on('connect', () => {
        console.log('[SocketIO] Connected');
        setIsConnected(true);
      });

      socketRef.current.on('connect_error', (error) => {
        console.error('[SocketIO] Connection error:', error);
        setIsConnected(false);
      });

      socketRef.current.on('disconnect', (reason) => {
        console.log('[SocketIO] Disconnected:', reason);
        setIsConnected(false);
      });

      socketRef.current.on('error', (error) => {
        console.error('[SocketIO] Socket error:', error);
        setIsConnected(false);
      });
    }

    return () => {
      if (socketRef.current) {
        socketRef.current.disconnect();
        socketRef.current = null;
      }
    };
  }, []);

  return (
    <SocketContext.Provider value={{ socket: socketRef.current, isConnected }}>
      {children}
    </SocketContext.Provider>
  );
}; 