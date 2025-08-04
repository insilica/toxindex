import React, { useState, useEffect } from 'react';

interface SessionTimeoutNotificationProps {
  isVisible: boolean;
  onExtendSession: () => void;
  onLogout: () => void;
  timeRemaining: number; // in seconds
}

export const SessionTimeoutNotification: React.FC<SessionTimeoutNotificationProps> = ({
  isVisible,
  onExtendSession,
  onLogout,
  timeRemaining
}) => {
  const [timeLeft, setTimeLeft] = useState(timeRemaining);

  useEffect(() => {
    if (isVisible && timeLeft > 0) {
      const timer = setInterval(() => {
        setTimeLeft(prev => {
          if (prev <= 1) {
            onLogout();
            return 0;
          }
          return prev - 1;
        });
      }, 1000);

      return () => clearInterval(timer);
    }
  }, [isVisible, timeLeft, onLogout]);

  useEffect(() => {
    setTimeLeft(timeRemaining);
  }, [timeRemaining]);

  if (!isVisible) return null;

  const minutes = Math.floor(timeLeft / 60);
  const seconds = timeLeft % 60;

  return (
    <div className="fixed top-4 right-4 z-50 bg-yellow-100 border border-yellow-400 text-yellow-800 px-4 py-3 rounded shadow-lg">
      <div className="flex items-center justify-between">
        <div className="flex-1">
          <h3 className="font-semibold">Session Timeout Warning</h3>
          <p className="text-sm">
            Your session will expire in {minutes}:{seconds.toString().padStart(2, '0')}
          </p>
        </div>
        <div className="flex space-x-2 ml-4">
          <button
            onClick={onExtendSession}
            className="bg-green-500 hover:bg-green-600 text-white px-3 py-1 rounded text-sm"
          >
            Extend Session
          </button>
          <button
            onClick={onLogout}
            className="bg-red-500 hover:bg-red-600 text-white px-3 py-1 rounded text-sm"
          >
            Logout Now
          </button>
        </div>
      </div>
    </div>
  );
}; 