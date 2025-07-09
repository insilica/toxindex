import React, { useEffect, useRef, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useEnvironment } from "../../context/EnvironmentContext";
import { IoServerOutline } from 'react-icons/io5';

interface EnvironmentSelectorProps {
  className?: string;
}

const EnvironmentSelector: React.FC<EnvironmentSelectorProps> = ({
  className = ''
}) => {
  const navigate = useNavigate();
  const selectRef = useRef<HTMLSelectElement>(null);
  const spanRef = useRef<HTMLSpanElement>(null);
  const [selectWidth, setSelectWidth] = useState(120);
  const [focused, setFocused] = useState(false);
  const { environments, selectedEnv, setSelectedEnv, loadingEnvironments } = useEnvironment();

  // Handle clicking outside to clear focus
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (selectRef.current && !selectRef.current.contains(event.target as Node)) {
        // Force blur and clear focus state
        selectRef.current.blur();
        setFocused(false);
        // Additional timeout to ensure state is cleared
        setTimeout(() => {
          setFocused(false);
        }, 10);
      }
    };

    document.addEventListener('click', handleClickOutside);
    return () => {
      document.removeEventListener('click', handleClickOutside);
    };
  }, []);

  // Update width based on selected content
  useEffect(() => {
    if (spanRef.current) {
      const measuredWidth = spanRef.current.offsetWidth;
      // Ensure minimum width of 200px, and add padding
      const calculatedWidth = Math.max(measuredWidth + 40, 300);
      setSelectWidth(calculatedWidth);
    }
  }, [selectedEnv, environments]);

  const handleChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const value = e.target.value;
    if (value === "__add__") {
      navigate("/settings/environments/create");
    } else if (value === "__manage__") {
      navigate("/settings/environments");
    } else {
      setSelectedEnv(value);
      navigate(`/environment/${value}`);
    }
  };

  return (
    <div className="relative group flex items-center gap-0" style={{ minWidth: 230, maxWidth: 400 }}>
      <IoServerOutline className="text-blue-400" style={{ fontSize: '1.6rem', flexShrink: 0 }} title="Environment" />
      {/* Hidden span to measure width */}
      <span
        ref={spanRef}
        style={{
          position: "absolute",
          visibility: "hidden",
          whiteSpace: "nowrap",
          fontWeight: "bold",
          fontSize: "0.875rem",
          fontFamily: "inherit",
          padding: "0 16px"
        }}
      >
        {selectedEnv === "__add__"
          ? "+ Add environment"
          : selectedEnv === "__manage__"
            ? "⚙ Manage environments"
            : (environments ?? []).find(e => e.environment_id === selectedEnv)?.title || 
              (environments && environments.length > 0 ? environments[0].title : "Select Environment")}
      </span>
      <select
        ref={selectRef}
        className={`font-bold text-white text-base leading-tight pl-2 pr-4 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-black group-hover:bg-opacity-60 group-hover:border group-hover:border-gray-700 group-hover:pl-4 group-hover:pr-12 group-hover:cursor-pointer ${focused ? 'bg-black bg-opacity-60 border border-gray-700 pl-4 pr-12' : ''} w-full ${loadingEnvironments ? 'opacity-50' : ''} ${className}`}
        style={{
          width: `${selectWidth}px`,
          minWidth: "100px",
          maxWidth: "300px",
          height: '2.5rem',
          borderRadius: '999px',
          position: 'relative',
          zIndex: 1,
          backgroundColor: 'transparent',
          cursor: 'pointer',
          paddingRight: '2.2rem',
          lineHeight: '1.5rem',
        }}
        value={selectedEnv || ''}
        onChange={handleChange}
        onFocus={() => setFocused(true)}
        disabled={loadingEnvironments}
      >
        {(environments ?? []).length > 0 && (environments ?? []).map(env => (
          <option key={env.environment_id} value={env.environment_id} style={{ paddingLeft: '1rem' }}>
            {env.title}
          </option>
        ))}
        <option value="__add__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>
          {loadingEnvironments ? '+ Add environment' : '+ Add environment'}
        </option>
        <option value="__manage__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>&#9881; Manage environments</option>
      </select>
      <span
        className="pointer-events-none absolute right-2 top-1/2 transform -translate-y-1/2 text-gray-400 group-hover:opacity-100 group-focus-within:opacity-100 opacity-0 transition-opacity duration-150"
        style={{ right: '1.2rem' }}
      >
        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
          <path d="M7 10l5 5 5-5" />
        </svg>
      </span>
    </div>
  );
};

export default EnvironmentSelector; 