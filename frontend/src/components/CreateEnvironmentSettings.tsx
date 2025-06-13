import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";

interface User {
  user_id: string;
  email: string;
}

interface CreateEnvironmentSettingsProps {
  refetchEnvironments: () => void;
}

const CreateEnvironmentSettings: React.FC<CreateEnvironmentSettingsProps> = ({ refetchEnvironments }) => {
  const [title, setTitle] = useState("");
  const [description, setDescription] = useState("");
  const [users, setUsers] = useState<User[]>([]);
  const [selectedUsers, setSelectedUsers] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    // Fetch all users for sharing
    fetch("/api/users", { credentials: "include" })
      .then(res => res.json())
      .then(data => setUsers(data.users || []));
  }, []);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setSuccess(null);
    if (!title.trim()) {
      setError("Environment name is required.");
      return;
    }
    const res = await fetch("/api/environment/new", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      credentials: "include",
      body: JSON.stringify({ title, description, sharing: selectedUsers }),
    });
    if (res.ok) {
      setSuccess("Environment created successfully!");
      refetchEnvironments();
      setTimeout(() => navigate("/settings/environments"), 1200);
    } else {
      setError("Failed to create environment.");
    }
  };

  const handleUserSelect = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const options = Array.from(e.target.selectedOptions, option => option.value);
    setSelectedUsers(options);
  };

  return (
    <div className="flex flex-col items-center min-h-screen min-w-full w-full h-screen pt-16" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-lg p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="font-bold text-white mb-6" style={{ fontSize: '1rem' }}>Create Environment</h1>
        <form onSubmit={handleSubmit} className="flex flex-col gap-4">
          <div>
            <label className="block text-white mb-1">Environment Name</label>
            <input
              type="text"
              value={title}
              onChange={e => setTitle(e.target.value)}
              className="w-full p-2 rounded border text-white"
              style={{ background: 'rgba(60,40,80,0.7)', borderColor: '#6b7280' }}
              placeholder="e.g. Lab workspace"
              required
            />
          </div>
          <div>
            <label className="block text-white mb-1">Description</label>
            <textarea
              value={description}
              onChange={e => setDescription(e.target.value)}
              className="w-full p-2 rounded border text-white"
              style={{ background: 'rgba(60,40,80,0.7)', borderColor: '#6b7280' }}
              rows={3}
              placeholder="Describe the purpose or scope..."
            />
          </div>
          <div>
            <label className="block text-white mb-1">Sharing (select users to share with)</label>
            <select
              multiple
              value={selectedUsers}
              onChange={handleUserSelect}
              className="w-full p-2 rounded border text-white"
              style={{ background: 'rgba(60,40,80,0.7)', borderColor: '#6b7280' }}
            >
              {users.map(user => (
                <option key={user.user_id} value={user.user_id}>{user.email}</option>
              ))}
            </select>
          </div>
          {error && <div className="text-red-500 text-sm">{error}</div>}
          {success && <div className="text-green-500 text-sm">{success}</div>}
          <button
            type="submit"
            className="mt-2 px-6 py-2 rounded-full font-semibold bg-green-600 text-white hover:bg-green-700 transition"
          >
            Create
          </button>
        </form>
      </div>
    </div>
  );
};

export default CreateEnvironmentSettings; 