import React from 'react';
import { FaSpinner } from 'react-icons/fa';

interface LoadingSpinnerProps {
  size?: 'small' | 'medium' | 'large';
  className?: string;
  text?: string;
}

const LoadingSpinner: React.FC<LoadingSpinnerProps> = ({
  size = 'medium',
  className = '',
  text = 'Loading...'
}) => {
  const sizeClass = {
    small: 'text-xl',
    medium: 'text-3xl',
    large: 'text-4xl'
  }[size];

  return (
    <div className={`flex items-center justify-center flex-col gap-3 ${className}`}>
      <FaSpinner className={`animate-spin ${sizeClass} text-gray-500`} />
      {text && <span className="text-gray-400">{text}</span>}
    </div>
  );
};

export default LoadingSpinner; 