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
  const [, setUsers] = useState<User[]>([]);
  const [selectedUsers, ] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState<string | null>(null);
  const [typedSharing, setTypedSharing] = useState("");
  const [typingSharing, setTypingSharing] = useState(true);
  const SHARING_PLACEHOLDER = "  Sharing coming soon... 🚧";
  const [typedName, setTypedName] = useState("");
  const [typingName, setTypingName] = useState(true);
  const NAME_PLACEHOLDER = "  e.g. Lab workspace 🧪";
  const [typedDesc, setTypedDesc] = useState("");
  const [typingDesc, setTypingDesc] = useState(true);
  const DESC_PLACEHOLDER = "  Describe the purpose or scope... 📝";
  const [nameFocused, setNameFocused] = useState(false);
  const [descFocused, setDescFocused] = useState(false);
  const navigate = useNavigate();

  useEffect(() => {
    // Fetch all users for sharing
    fetch("/api/users", { credentials: "include" })
      .then(res => res.json())
      .then(data => setUsers(data.users || []));
  }, []);

  useEffect(() => {
    let i = 0;
    setTypedSharing("");
    setTypingSharing(true);
    const interval = setInterval(() => {
      setTypedSharing(SHARING_PLACEHOLDER.slice(0, i + 1));
      i++;
      if (i === SHARING_PLACEHOLDER.length) {
        clearInterval(interval);
        setTypingSharing(false);
      }
    }, 60);
    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    let i = 0;
    setTypedName("");
    setTypingName(true);
    const interval = setInterval(() => {
      setTypedName(NAME_PLACEHOLDER.slice(0, i + 1));
      i++;
      if (i === NAME_PLACEHOLDER.length) {
        clearInterval(interval);
        setTypingName(false);
      }
    }, 50);
    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    let i = 0;
    setTypedDesc("");
    setTypingDesc(true);
    const interval = setInterval(() => {
      setTypedDesc(DESC_PLACEHOLDER.slice(0, i + 1));
      i++;
      if (i === DESC_PLACEHOLDER.length) {
        clearInterval(interval);
        setTypingDesc(false);
      }
    }, 40);
    return () => clearInterval(interval);
  }, []);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setSuccess(null);
    if (!title.trim()) {
      setError("Environment name is required.");
      return;
    }
    const res = await fetch("/api/environments", {
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

  // const handleUserSelect = (e: React.ChangeEvent<HTMLSelectElement>) => {
  //   const options = Array.from(e.target.selectedOptions, option => option.value);
  //   setSelectedUsers(options);
  // };

  return (
    <div className="flex flex-col items-center min-h-screen min-w-full w-full h-screen pt-52" 
    style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-lg p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="font-bold text-white mb-6" style={{ fontSize: '1rem' }}>Create Environment</h1>
        <form onSubmit={handleSubmit} className="flex flex-col gap-4">
          <div>
            <label className="block text-white mb-1">Environment Name</label>
            <div style={{ position: 'relative' }}>
              <input
                type="text"
                value={title}
                onChange={e => setTitle(e.target.value)}
                onFocus={() => setNameFocused(true)}
                onBlur={() => setNameFocused(false)}
                className="w-full p-2 rounded border text-white bg-transparent"
                style={{ background: 'rgba(60,40,80,0.7)', borderColor: '#6b7280', position: 'relative', zIndex: 2 }}
                required
              />
              {(!title && !nameFocused) && (
                <span style={{
                  position: 'absolute',
                  left: 12,
                  top: '50%',
                  transform: 'translateY(-50%)',
                  color: '#a3a3a3',
                  fontFamily: 'monospace',
                  fontSize: '1.15em',
                  letterSpacing: '0.01em',
                  pointerEvents: 'none',
                  zIndex: 1,
                  opacity: 0.85
                }}>
                  {typedName}
                  {typingName && <span style={{ color: '#a3a3a3', fontWeight: 'bold', animation: 'blink 1s steps(1) infinite' }}>|</span>}
                </span>
              )}
            </div>
          </div>
          <div>
            <label className="block text-white mb-1">Description</label>
            <div style={{ position: 'relative' }}>
              <textarea
                value={description}
                onChange={e => setDescription(e.target.value)}
                onFocus={() => setDescFocused(true)}
                onBlur={() => setDescFocused(false)}
                className="w-full p-2 rounded border text-white bg-transparent"
                style={{ background: 'rgba(60,40,80,0.7)', borderColor: '#6b7280', position: 'relative', zIndex: 2 }}
                rows={3}
              />
              {(!description && !descFocused) && (
                <span style={{
                  position: 'absolute',
                  left: 12,
                  top: 12,
                  color: '#a3a3a3',
                  fontFamily: 'monospace',
                  fontSize: '1.15em',
                  letterSpacing: '0.01em',
                  pointerEvents: 'none',
                  zIndex: 1,
                  opacity: 0.85
                }}>
                  {typedDesc}
                  {typingDesc && <span style={{ color: '#a3a3a3', fontWeight: 'bold', animation: 'blink 1s steps(1) infinite' }}>|</span>}
                </span>
              )}
            </div>
          </div>
          <div>
            <label className="block text-white mb-1">Sharing (select users to share with)</label>
            <div className="w-full p-2 rounded border text-gray-400 bg-gray-800 opacity-60 cursor-not-allowed relative" style={{ borderColor: '#6b7280', minHeight: 44, display: 'flex', alignItems: 'center' }}>
              <span style={{ color: '#a3a3a3', fontFamily: 'monospace', fontSize: '1.15em', letterSpacing: '0.01em', minHeight: 24, opacity: 0.85 }}>
                {typedSharing}
                {typingSharing && <span style={{ color: '#a3a3a3', fontWeight: 'bold', animation: 'blink 1s steps(1) infinite' }}>|</span>}
              </span>
            </div>
            <div className="text-gray-400 text-xs mt-2" style={{ fontStyle: 'italic' }}>🔧 Sharing is under construction. Stay tuned for collaboration features!</div>
            <style>{`
              @keyframes blink {
                0%, 100% { opacity: 1; }
                50% { opacity: 0; }
              }
            `}</style>
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