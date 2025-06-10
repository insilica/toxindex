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
    <div className="max-w-3xl mx-auto mt-12 p-6 bg-white rounded shadow">
      <h1 className="text-2xl font-bold mb-6">Environments</h1>
      <table className="w-full border-collapse">
        <thead>
          <tr className="bg-gray-200">
            <th className="p-2 text-left">Name</th>
            <th className="p-2 text-left">Number of Tasks</th>
            <th className="p-2 text-left">Creator</th>
            <th className="p-2 text-left">Created At</th>
          </tr>
        </thead>
        <tbody>
          {environments.map(env => (
            <tr key={env.environment_id} className="border-b">
              <td className="p-2">{env.title}</td>
              <td className="p-2">0</td>
              <td className="p-2">{env.user_id}</td>
              <td className="p-2">{new Date(env.created_at).toLocaleDateString(undefined, { month: '2-digit', day: '2-digit', year: 'numeric' })}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default SettingsEnvironments; 