import React, { useState } from "react";
import { useNavigate } from "react-router-dom";

const CreateEnvironment: React.FC = () => {
  const [title, setTitle] = useState("");
  const [description, setDescription] = useState("");
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);

    const res = await fetch("/api/environment/new", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      credentials: "include",
      body: JSON.stringify({ title, description }),
    });

    if (res.ok) {
      const data = await res.json();
      if (data.environment_id) {
        navigate(`/environment/${data.environment_id}`);
      }
    } else {
      setError("Failed to create environment.");
    }
  };

  return (
    <form onSubmit={handleSubmit} className="max-w-md mx-auto p-6 bg-white rounded shadow">
      <h2 className="text-xl font-bold mb-4">Create Environment</h2>
      <input
        type="text"
        placeholder="Environment name"
        value={title}
        onChange={e => setTitle(e.target.value)}
        required
        className="w-full mb-3 p-2 border rounded"
      />
      <textarea
        placeholder="Description (optional)"
        value={description}
        onChange={e => setDescription(e.target.value)}
        className="w-full mb-3 p-2 border rounded"
      />
      {error && <div className="text-red-500 mb-2">{error}</div>}
      <button
        type="submit"
        className="w-full bg-blue-600 text-white py-2 rounded hover:bg-blue-700"
      >
        Create
      </button>
    </form>
  );
};

export default CreateEnvironment; 