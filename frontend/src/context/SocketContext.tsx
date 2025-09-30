import React, { createContext, useContext, useEffect, useState } from 'react';
import { io, Socket } from 'socket.io-client';

interface SocketContextType {
  socket: Socket | null;
  isConnected: boolean;
  connect: () => void;
  disconnect: () => void;
  isConnecting: boolean;
}

const SocketContext = createContext<SocketContextType>({
  socket: null,
  isConnected: false,
  connect: () => {},
  disconnect: () => {},
  isConnecting: false,
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
  console.log('[SocketProvider] SocketProvider mounting');
  const [socket, setSocket] = useState<Socket | null>(null);
  const [isConnected, setIsConnected] = useState(false);
  const [isConnecting, setIsConnecting] = useState(false);

  useEffect(() => {
    console.log('[SocketProvider] SocketProvider mounted');
    return () => {
      console.log('[SocketProvider] SocketProvider unmounting');
    };
  }, []);

  const connect = () => {
    // Don't connect if already connecting or connected
    if (isConnecting || isConnected || socket) {
      console.log('[SocketIO] Connection already in progress or established');
      return;
    }

    // Check authentication first
    fetch("/api/users/me", { credentials: "include", cache: "no-store" })
      .then(res => {
        if (res.status === 200) {
          console.log('[SocketIO] User authenticated, connecting...');
          establishConnection();
        } else {
          console.log('[SocketIO] User not authenticated');
        }
      })
      .catch(() => {
        console.log('[SocketIO] Authentication check failed');
      });
  };

  const establishConnection = () => {
    if (socket) {
      console.log('[SocketIO] Socket already exists, disconnecting first');
      socket.disconnect();
      setSocket(null);
    }

    console.log('[SocketIO] Creating new socket connection...');
    setIsConnecting(true);
    const s = io('/', {
      withCredentials: true,
      transports: ['websocket', 'polling'],
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
    setSocket(s);

    // Add connection event handlers
    s.on('connect', () => {
      console.log('[SocketIO] Connected');
      setIsConnected(true);
      setIsConnecting(false);
    });

    s.on('connect_error', (error) => {
      console.error('[SocketIO] Connection error:', error);
      setIsConnected(false);
      setIsConnecting(false);
    });

    s.on('disconnect', (reason) => {
      console.log('[SocketIO] Disconnected:', reason);
      setIsConnected(false);
      setIsConnecting(false);
    });

    s.on('error', (error) => {
      console.error('[SocketIO] Socket error:', error);
      setIsConnected(false);
      setIsConnecting(false);
    });
  };

  const disconnect = () => {
    if (socket) {
      console.log('[SocketIO] Disconnecting socket...');
      socket.disconnect();
      setSocket(null);
      setIsConnected(false);
      setIsConnecting(false);
    }
  };

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      disconnect();
    };
  }, []);

  return (
    <SocketContext.Provider value={{ 
      socket, 
      isConnected, 
      connect, 
      disconnect,
      isConnecting
    }}>
      {children}
    </SocketContext.Provider>
  );
}; 