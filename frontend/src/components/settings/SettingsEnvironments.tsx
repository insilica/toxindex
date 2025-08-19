import React, {useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { useEnvironment } from "../../context/EnvironmentContext";
import HomeButton from '../shared/HomeButton';

interface Environment {
  environment_id: string;
  title: string;
  user_id: string;
  created_at: string;
  task_count?: number;
  file_count?: number;
}

const SettingsEnvironments: React.FC = () => {
  const [modalOpen, setModalOpen] = useState(false);
  const [envToDelete, setEnvToDelete] = useState<Environment | null>(null);
  const navigate = useNavigate();
  const { environments, setEnvironments, refetchEnvironments } = useEnvironment();

  useEffect(() => {
    // Fetch environments when component mounts
    refetchEnvironments();
  }, []); // Only run once on mount

  const handleDelete = (envId: string) => {
    fetch(`/api/environments/${envId}`, {
      method: 'DELETE',
      credentials: 'include',
      
    })
      .then(res => res.json())
      .then(data => {
        if (data.success) {
          if (setEnvironments) {
            setEnvironments(envs => envs.filter(env => env.environment_id !== envId));
          }
          refetchEnvironments();
        } else {
          alert('Failed to delete environment.');
        }
      })
      .catch(() => alert('Failed to delete environment.'));
    setModalOpen(false);
    setEnvToDelete(null);
  };

  return (
    <div className="max-w-6xl mx-auto p-6 relative" style={{ paddingLeft: '8rem' }}>
      <div className="w-full max-w-3xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        {/* Title and Add Button Row */}
        <div className="flex items-center justify-between mb-4">
          <h1 className="font-bold text-white" style={{marginBottom: 0, fontSize: '1.2rem'}}>Environments</h1>
        </div>
        {/* Underline */}
        <div className="w-full border-b border-gray-400 mb-4"></div>
        {/* Add Button Row */}
        <div className="flex justify-end items-center mb-8">
          <button
            type="button"
            onClick={() => navigate('/settings/environments/create')}
            className="px-6 font-semibold transition"
            style={{
              fontWeight: 600,
              fontSize: '.8rem',
              paddingTop: 0,
              paddingBottom: 0,
              lineHeight: 1.1,
              boxShadow: '0 1px 4px 0 rgba(0,0,0,0.04)',
              background: '#f3f4f6',
              color: '#111',
              borderRadius: '9999px',
            }}
            onMouseOver={e => (e.currentTarget.style.background = '#e5e7eb')}
            onMouseOut={e => (e.currentTarget.style.background = '#f3f4f6')}
          >
            <span style={{ fontSize: '2.5em', fontWeight: 200, verticalAlign: 'middle' }}>+  </span> 
            <span style={{ fontSize: '.9rem', fontWeight: 600, verticalAlign: 'middle' }}>Create environment</span>
          </button>
        </div>
        <table className="w-full border-collapse text-white">
          <thead>
            <tr className="bg-transparent">
              <th className="p-2 text-left">Name</th>
              <th className="p-2 text-left">Number of tasks</th>
              <th className="p-2 text-left">Number of files</th>
              <th className="p-2 text-left">Creator</th>
              <th className="p-2 text-left">Created at</th>
            </tr>
          </thead>
          <tbody>
            {environments.map(env => (
              <tr
                key={env.environment_id}
                className={`border-b border-gray-700 transition-colors hover:bg-purple-900/30 cursor-pointer`}
                onClick={() => navigate(`/environments/details?env=${env.environment_id}`)}
              >
                <td className="py-3 px-4 text-white text-lg font-medium">
                  {env.title}
                </td>
                <td className="py-3 px-4 text-white text-lg font-medium">{env.task_count ?? 0}</td>
                <td className="py-3 px-4 text-white text-lg font-medium">{env.file_count ?? 0}</td>
                <td className="py-3 px-4 text-white text-lg font-medium">{env.user_id}</td>
                <td className="py-3 px-4 text-gray-400 text-sm">
                  {env.created_at && new Date(env.created_at).toLocaleDateString(undefined, { month: 'short', day: 'numeric', year: 'numeric' })}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {/* Modal overlay for delete confirmation */}
        {modalOpen && envToDelete && (
          <div
            style={{
              position: 'fixed',
              top: 0,
              left: 0,
              width: '100vw',
              height: '100vh',
              background: 'rgba(0,0,0,0.55)',
              zIndex: 1000,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
            }}
            onClick={() => { setModalOpen(false); setEnvToDelete(null); }}
          >
            <div
              style={{
                background: '#181c1f',
                borderRadius: 16,
                padding: '2.5rem 2.5rem 2rem 2.5rem',
                minWidth: 320,
                boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)',
                border: '1.5px solid #b91c1c',
                position: 'relative',
              }}
              onClick={e => e.stopPropagation()}
            >
              <div className="text-lg text-white font-semibold mb-4" style={{ textAlign: 'center' }}>
                Are you sure you want to delete <span className="text-red-400">{envToDelete.title}</span>?
              </div>
              <div className="flex justify-center gap-6 mt-6">
                <button
                  className="px-6 py-2 rounded bg-gray-700 text-white hover:bg-gray-800 font-semibold"
                  onClick={() => { setModalOpen(false); setEnvToDelete(null); }}
                >
                  Cancel
                </button>
                <button
                  className="px-6 py-2 rounded bg-red-700 text-white hover:bg-red-800 font-semibold"
                  onClick={() => handleDelete(envToDelete.environment_id)}
                >
                  Delete
                </button>
              </div>
            </div>
          </div>
        )}
      </div>
      
      {/* Home button positioned relative to settings content */}
      <div 
        className="absolute transition-all duration-300 z-50"
        style={{
          left: '2rem',
          top: '1.5rem',
          border: 'none',
          padding: 0
        }}
      >
        <HomeButton
          color="#16a34a"
          hoverColor="#2563eb"
          aria-label="Go back"
        />
      </div>
    </div>
  );
};

export default SettingsEnvironments; 