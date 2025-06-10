import React, { useEffect, useState } from "react";

interface Environment {
  environment_id: string;
  title: string;
  user_id: string;
  created_at: string;
}

const SettingsEnvironments: React.FC = () => {
  const [environments, setEnvironments] = useState<Environment[]>([]);

  useEffect(() => {
    fetch("/api/environments", { credentials: "include" })
      .then(res => res.json())
      .then(data => setEnvironments(data.environments || []));
  }, []);

  return (
    <div className="flex items-center justify-center min-h-screen min-w-full w-full h-screen" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-3xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-3xl font-bold mb-8 text-white text-center">Environments</h1>
        <table className="w-full border-collapse text-white">
          <thead>
            <tr className="bg-transparent">
              <th className="p-2 text-left">Name</th>
              <th className="p-2 text-left">Number of Tasks</th>
              <th className="p-2 text-left">Creator</th>
              <th className="p-2 text-left">Created At</th>
            </tr>
          </thead>
          <tbody>
            {environments.map(env => (
              <tr key={env.environment_id} className="border-b border-gray-700">
                <td className="p-2">{env.title}</td>
                <td className="p-2">0</td>
                <td className="p-2">{env.user_id}</td>
                <td className="p-2">{new Date(env.created_at).toLocaleDateString(undefined, { month: '2-digit', day: '2-digit', year: 'numeric' })}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default SettingsEnvironments; 