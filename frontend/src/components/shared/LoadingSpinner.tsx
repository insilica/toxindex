import React from 'react';
import { FaSpinner } from 'react-icons/fa';

interface LoadingSpinnerProps {
  size?: 'small' | 'medium' | 'large';
  className?: string;
  text?: string;
  showTimer?: boolean;
  startTime?: number; // UNIX timestamp in seconds or ms
}

const LoadingSpinner: React.FC<LoadingSpinnerProps> = ({
  size = 'medium',
  className = '',
  text = '',
  showTimer = false,
  startTime
}) => {
  const sizeClass = {
    small: 'text-xl',
    medium: 'text-3xl',
    large: 'text-4xl'
  }[size];

  // If startTime is provided, calculate elapsed seconds from it
  const getInitialSeconds = () => {
    if (!startTime) return 0;
    const now = Date.now();
    // If startTime is in seconds, convert to ms
    const startMs = startTime > 1e12 ? startTime : startTime * 1000;
    return Math.floor((now - startMs) / 1000);
  };

  const [seconds, setSeconds] = React.useState(getInitialSeconds());
  React.useEffect(() => {
    if (!showTimer) return;
    const update = () => setSeconds(getInitialSeconds());
    update(); // set immediately
    const interval = setInterval(update, 1000);
    return () => clearInterval(interval);
  }, [showTimer, startTime]);

  return (
    <div className={`flex items-center justify-center flex-col gap-3 ${className}`}>
      <div className="flex items-center gap-2">
        <FaSpinner className={`animate-spin ${sizeClass} text-gray-500`} />
        {showTimer && (
          <span className="text-xs text-gray-400 font-mono" aria-label={`It's been ${seconds} seconds`}>
            {seconds}secs
          </span>
        )}
      </div>
      {text && <span className="text-gray-400">{text}</span>}
    </div>
  );
};

export default LoadingSpinner; 