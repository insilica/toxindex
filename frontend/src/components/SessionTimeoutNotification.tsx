import React, { useState, useEffect, useRef } from 'react';

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
  const [localTimeRemaining, setLocalTimeRemaining] = useState(timeRemaining);
  const countdownIntervalRef = useRef<number | null>(null);

  // Update local time when server time changes
  useEffect(() => {
    setLocalTimeRemaining(timeRemaining);
  }, [timeRemaining]);

  // Set up local countdown timer
  useEffect(() => {
    if (countdownIntervalRef.current) {
      clearInterval(countdownIntervalRef.current);
    }

    if (isVisible && localTimeRemaining > 0) {
      countdownIntervalRef.current = window.setInterval(() => {
        setLocalTimeRemaining(prev => {
          const newValue = Math.max(0, prev - 1);
          return newValue;
        });
      }, 1000);
    }

    return () => {
      if (countdownIntervalRef.current) {
        clearInterval(countdownIntervalRef.current);
      }
    };
  }, [isVisible, localTimeRemaining]);

  // Auto-logout when time reaches 0
  useEffect(() => {
    if (isVisible && localTimeRemaining <= 0) {
      onLogout();
    }
  }, [isVisible, localTimeRemaining, onLogout]);

  if (!isVisible) return null;

  const minutes = Math.floor(localTimeRemaining / 60);
  const seconds = localTimeRemaining % 60;

  return (
    <>
      {/* Backdrop with foggy effect */}
      <div className="fixed inset-0 backdrop-blur-sm z-40" />
      
      {/* Modal in center */}
      <div className="fixed inset-0 z-50 flex items-center justify-center">
        <div 
          className="bg-white border border-gray-300 text-gray-800 px-6 py-4 rounded-lg shadow-xl max-w-sm mx-4 cursor-pointer hover:bg-gray-50 transition-colors"
          onClick={onExtendSession}
        >
          <div className="text-center">
            <h3 className="text-lg font-semibold mb-2">Session Timeout</h3>
            <p className="text-sm text-gray-600 mb-4">
              Your session will expire in {minutes}:{seconds.toString().padStart(2, '0')}
            </p>
            <p className="text-sm text-gray-500">
              Click anywhere to extend your session
            </p>
          </div>
        </div>
      </div>
    </>
  );
}; 