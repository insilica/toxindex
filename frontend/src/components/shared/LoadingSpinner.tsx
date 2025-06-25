import React from 'react';
import { FaSpinner } from 'react-icons/fa';

interface LoadingSpinnerProps {
  size?: 'small' | 'medium' | 'large';
  className?: string;
  text?: string;
  showTimer?: boolean;
}

const LoadingSpinner: React.FC<LoadingSpinnerProps> = ({
  size = 'medium',
  className = '',
  text = '',
  showTimer = false
}) => {
  const sizeClass = {
    small: 'text-xl',
    medium: 'text-3xl',
    large: 'text-4xl'
  }[size];

  const [seconds, setSeconds] = React.useState(0);
  React.useEffect(() => {
    if (!showTimer) return;
    const interval = setInterval(() => setSeconds(s => s + 1), 1000);
    return () => clearInterval(interval);
  }, [showTimer]);

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