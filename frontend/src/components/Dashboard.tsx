import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";

interface Environment {
  environment_id: string;
  title: string;
}

const Dashboard: React.FC = () => {
  const [environments, setEnvironments] = useState<Environment[]>([]);
  const [selectedEnv, setSelectedEnv] = useState<string>("");
  const [showAdd, setShowAdd] = useState(false);
  const [newEnvTitle, setNewEnvTitle] = useState("");
  const [chatInput, setChatInput] = useState("is aspirin hepatotoxic?");
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    fetch("/environments", { credentials: "include" })
      .then(res => res.json())
      .then(data => {
        setEnvironments(data.environments || []);
        if (data.environments && data.environments.length > 0) {
          setSelectedEnv(data.environments[0].environment_id);
        }
      });
  }, []);

  const handleAddEnvironment = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    const res = await fetch("/api/environment/new", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      credentials: "include",
      body: JSON.stringify({ title: newEnvTitle }),
    });
    if (res.ok) {
      const data = await res.json();
      setShowAdd(false);
      setNewEnvTitle("");
      // Refresh environments
      fetch("/environments", { credentials: "include" })
        .then(res => res.json())
        .then(data => {
          setEnvironments(data.environments || []);
          if (data.environments && data.environments.length > 0) {
            setSelectedEnv(data.environments[data.environments.length - 1].environment_id);
          }
        });
    } else {
      setError("Failed to create environment.");
    }
  };

  const handleChatSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    // You can handle chat input here (send to backend, etc.)
    setChatInput("");
  };

  return (
    <div className="max-w-2xl mx-auto min-h-screen flex flex-col justify-center" style={{ background: '#10291a' }}>
      <h2 className="text-xl font-semibold mb-4 text-white">Toxindex Dashboard</h2>
      <form onSubmit={handleChatSubmit} className="relative flex flex-col gap-4 items-stretch">
        {/* Chat-style textarea with dropdown inside */}
        <div className="relative">
          <textarea
            rows={4}
            placeholder="Type your message..."
            value={chatInput}
            onChange={e => setChatInput(e.target.value)}
            className="w-full p-5 pr-40 text-lg rounded-2xl border border-gray-300 bg-white bg-opacity-80 text-black resize-none min-h-[96px] shadow-lg focus:outline-none focus:ring-2 focus:ring-green-400"
            style={{ minHeight: 120, fontFamily: 'inherit' }}
          />
          {/* Dropdown inside the textbox, top-right */}
          <div className="absolute top-3 right-3 z-10">
            <select
              className="p-1 px-2 rounded bg-white bg-opacity-60 text-black text-sm border border-gray-300 shadow-sm"
              style={{ minWidth: 120 }}
              value={selectedEnv}
              onChange={e => setSelectedEnv(e.target.value)}
            >
              {environments.map(env => (
                <option key={env.environment_id} value={env.environment_id}>
                  {env.title}
                </option>
              ))}
              <option value="__add__">+ Add environment</option>
            </select>
          </div>
        </div>
        {/* Add environment inline form */}
        {selectedEnv === "__add__" && (
          <form onSubmit={handleAddEnvironment} className="flex gap-2 mt-2">
            <input
              type="text"
              placeholder="New environment name"
              value={newEnvTitle}
              onChange={e => setNewEnvTitle(e.target.value)}
              required
              className="flex-1 p-2 border rounded bg-white bg-opacity-70 text-black"
            />
            <button type="submit" className="bg-blue-600 text-white px-4 py-2 rounded hover:bg-blue-700">
              Add
            </button>
            <button type="button" onClick={() => { setShowAdd(false); setSelectedEnv(environments[0]?.environment_id || ""); }} className="px-4 py-2 rounded border">
              Cancel
            </button>
          </form>
        )}
        <button type="submit" className="self-end bg-green-600 text-white px-6 py-3 rounded-2xl text-lg font-semibold hover:bg-green-700 mt-2 shadow">
          Send
        </button>
        {error && <div className="text-red-500 mt-2">{error}</div>}
      </form>
    </div>
  );
};

export default Dashboard; 