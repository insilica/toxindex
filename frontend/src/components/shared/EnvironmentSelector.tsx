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
  const { environments, selectedEnv, setSelectedEnv, loadingEnvironments } = useEnvironment();

  // Update width based on selected content
  useEffect(() => {
    if (spanRef.current) {
      setSelectWidth(spanRef.current.offsetWidth + 40); // add padding
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
    <div className="relative group flex items-center gap-2" style={{ minWidth: 180, maxWidth: 210 }}>
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
            : (environments ?? []).find(e => e.environment_id === selectedEnv)?.title || ""}
      </span>
      <select
        ref={selectRef}
        className={`font-bold text-white text-base leading-tight px-4 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-black group-hover:bg-opacity-60 group-hover:border group-hover:border-gray-700 group-hover:px-6 group-hover:pr-12 group-hover:cursor-pointer focus:bg-black focus:bg-opacity-60 focus:border focus:border-gray-700 focus:px-6 focus:pr-12 w-full ${loadingEnvironments ? 'opacity-50' : ''} ${className}`}
        style={{
          width: `${selectWidth}px`,
          minWidth: "50px",
          maxWidth: "260px",
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